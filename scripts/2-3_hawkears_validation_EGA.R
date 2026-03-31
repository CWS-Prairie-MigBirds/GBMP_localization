#Establishing Species-Specific Thresholds for CNN (HawkEars)
#By: Megan Edgar (Jan 2026)
#Modified: Erica Alex (Jan 2026)

#Adapted from the Methods of Tseng et al (2025)
#Setting BirdNET confidence thresholds: species-specific vs. universal approaches. J Ornithol 166, 1123–1135 (2025). 
#https://doi.org/10.1007/s10336-025-02260-w

#Libraries --------------------------------------------------------------------------------------------------------------------------------------------------------

library(dplyr)
library(purrr)
library(tidyr)
library(tuneR)
library(soundecology)
library(readr)
library(ggplot2)

#Load data --------------------------------------------------------------------------------------------------------------------------------------------------------

validation <- read.csv("H:/Localization/grsp_valid_0204.csv")  #validated dataset
labels <- read.csv("H:/Localization/hawkears output/HawkEars_labels.csv")     #full HawkEars run
audio_root <- "H:/Localization/raw recordings/GSA/GSA-trimmed/20250618_0803"                           #location of audio files with matching names

# validation2 <- read.csv("C:/Users/AlexE/cclo_trim/cclo_valid_0202.csv")
# 
# validation = validation1 %>%
#   bind_rows(validation2) %>%
#   distinct(filename, start_time, end_time, .keep_all=TRUE)

# validation2 <- read.csv("C:/Users/AlexE/grsp_trim/grsp_valid_xtra.csv")
# 
# validation = validation1 %>%
#   bind_rows(validation2) %>%
#   distinct(filename, start_time, end_time, .keep_all=TRUE)


#Tseng method - no acoustic indices ------------------------------------------------------------------------------------------------------------------------------

#set bins and target precision
min_conf <- 0.10
max_conf <- 1.00
step     <- 0.05
target_precision <- 0.90
fallback_threshold <- 0.95  # Tseng fallback when 0.9 can't be achieved 
threshold_grid <- seq(min_conf, fallback_threshold, by = step)

#clean validation data
val <- as_tibble(validation) %>%
  mutate(
    species = class_code,
    score   = as.numeric(score),
    tp = case_when(
      tolower(label) %in% c("yes","y","true","1","present") ~ 1L,
      tolower(label) %in% c("no","n","false","0","absent")  ~ 0L,
      TRUE ~ NA_integer_
    )
  ) %>%
  filter(!is.na(tp), !is.na(score), !is.na(species)) %>%
  filter(score >= min_conf, score <= max_conf)

#perform regression on true postives and HE score
val_models <- val %>%
  group_by(species) %>%
  nest() %>% #collapse each group into a list-column
  mutate(
    model = map(data, ~ glm(tp ~ score, family = binomial(), data = .x)), #logistic regression looped over each nested group
    n_val = map_int(data, nrow) #count number of validation obs for each species
  ) %>%
  select(species, model, n_val)

fits <- val_models

#Bin full labels into 0.05 confidence classes and count N_i per species (#detections in each bin)
lab <- as_tibble(labels) %>%
  mutate(
    species = class_code,       # must match val species choice
    score   = as.numeric(score)
  ) %>%
  filter(!is.na(score), !is.na(species)) %>%
  filter(score >= min_conf, score <= max_conf) %>%
  semi_join(fits %>% select(species), by = "species")  # only species  validated

# Fast numeric binning to the LOWER edge of the 0.05 bin:
lab_bins <- lab %>%
  mutate(
    bin_lower = min_conf + floor((score - min_conf) / step) * step,
    bin_lower = pmin(fallback_threshold, pmax(min_conf, bin_lower)),
    bin_mid   = bin_lower + step/2
  ) %>%
  count(species, bin_lower, bin_mid, name = "N_i")

#For each species/bin, predict TPR_i from GLM and compute NTP_i / NFP_i
#Then compute precision(T) over thresholds and pick minimum T with precision>=0.9

#function to compute threshold
#create table
compute_tseng_threshold <- function(df_bins, model) {
  df_bins <- df_bins %>%
    mutate(
      TPR_i = predict(model, newdata = data.frame(score = bin_mid), type = "response"),         #predict prob of true positive at score
      NTP_i = TPR_i * N_i,                                                                      #expected number of true positives
      NFP_i = (1 - TPR_i) * N_i                                                                 #expected number of false positives
    )
  
  curve <- purrr::map_dfr(threshold_grid, function(T) {
    kept <- df_bins %>% dplyr::filter(bin_lower >= T)                                           #keep bins above the threshold
    denom <- sum(kept$NTP_i + kept$NFP_i)                                                       #precision = expected TP/expected(TP+FP)
    prec  <- if (denom == 0) NA_real_ else sum(kept$NTP_i) / denom                              #handle for /0
    retained <- if (sum(df_bins$N_i) == 0) NA_real_ else sum(kept$N_i) / sum(df_bins$N_i)       #proportion of obs that are kept based on the threshold
    
    tibble::tibble(
      threshold = T,
      precision = prec,
      prop_retained = retained,
      n_retained = sum(kept$N_i)
    )
  })
  #select lowest threshold that meets target precision  
  hit <- curve %>%
    dplyr::filter(!is.na(precision), precision >= target_precision) %>%
    dplyr::arrange(threshold) %>%
    dplyr::slice(1)
  
  status <- "hit_target"
  if (nrow(hit) == 0) {
    hit <- curve %>% dplyr::filter(threshold == fallback_threshold)
    status <- "fallback"
  }
  
  hit <- hit %>% dplyr::mutate(status = status)
  
  list(summary = hit, curve = curve, bins = df_bins)
}
#join fitted models to binned validation data, group by species, compute threshold for each species
tables <- fits %>%
  dplyr::select(species, model, n_val) %>%
  dplyr::left_join(lab_bins %>% dplyr::group_by(species) %>%
                     dplyr::summarise(total_detections = sum(N_i), .groups = "drop"),
                   by = "species") %>%
  dplyr::left_join(lab_bins, by = "species") %>%
  dplyr::group_by(species, model, n_val, total_detections) %>%
  dplyr::summarise(
    out = list(
      compute_tseng_threshold(
        dplyr::cur_data_all() %>% dplyr::select(bin_lower, bin_mid, N_i),
        model[[1]]
      )
    ),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    chosen_threshold = purrr::map_dbl(out, ~ .x$summary$threshold),
    precision_at_threshold = purrr::map_dbl(out, ~ .x$summary$precision),
    prop_retained = purrr::map_dbl(out, ~ .x$summary$prop_retained),
    n_retained = purrr::map_int(out, ~ .x$summary$n_retained),
    status = purrr::map_chr(out, ~ .x$summary$status),
    n_dropped = total_detections - n_retained
  )

threshold_table = tables %>%
  dplyr::select(
    species,
    n_val,
    total_detections,
    chosen_threshold,
    precision_at_threshold,
    prop_retained,
    n_retained,
    n_dropped,
    status
  ) %>%
  dplyr::arrange(dplyr::desc(chosen_threshold), species)

threshold_table #one row per species

curves_table <- tables %>%
  dplyr::select(species, out) %>%
  tidyr::unnest_wider(out) %>%
  dplyr::select(species, curve) %>%
  tidyr::unnest(curve)

curves_table #many rows per species x thresholds

#Plot species curve with ggplot
plot_species_gg_onepanel <- function(species_code, target = 0.90) {
  
  d <- curves_table %>%
    dplyr::filter(species == species_code) %>%
    dplyr::arrange(threshold)
  
  if(nrow(d)==0) stop("Species curve not found")
  
  thr_star <- threshold_table %>%
    dplyr::filter(species == species_code) %>%
    dplyr::pull(chosen_threshold)
  
  d_long <- d %>%
    dplyr::select(threshold, precision, prop_retained) %>%
    tidyr:: pivot_longer(
      cols = c(precision, prop_retained),
      names_to = "metric",
      values_to = "value"
    ) %>%
    dplyr::mutate(
      metric = dplyr::recode(
        metric,
        precision = "Precision",
        prop_retained = 'Proportion retained'
      )
    )
  
  ggplot2::ggplot(d_long, ggplot2::aes(x = threshold, y = value, linetype = metric)) +
    ggplot2::geom_line(linewidth = 1.1) +
    ggplot2::geom_point(size = 1.6, alpha = 0.85) +
    ggplot2::geom_hline(yintercept = target, linetype = "dashed", linewidth = 0.8) +
    ggplot2::geom_vline(xintercept = thr_star, linetype = "dotdash", linewidth = 0.8) +
    ggplot2::scale_y_continuous(
      limits = c(0, 1),
      breaks = seq(0, 1, by = 0.2),
      name = "Precision of retained detections",
      sec.axis = ggplot2::sec_axis(~ ., name = "Proportion of detections retained")
    ) +
    ggplot2::scale_x_continuous(name = "Confidence threshold (T)") +
    ggplot2::scale_linetype_manual(values = c("Precision" = "solid", "Proportion retained" = "solid")) +
    ggplot2::labs(
      title = paste0(species_code, " — Tseng-style threshold (precision + retained)"),
      subtitle = paste0("Chosen threshold T* = ", format(thr_star, digits = 3),
                        " | Target precision = ", target),
      linetype = NULL
    ) +
    ggplot2::annotate(
      "text",
      x = thr_star,
      y = 0.98,
      label = paste0("T* = ", format(thr_star, digits = 3)),
      hjust = -0.05,
      vjust = 1,
      size = 3.6
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold"),
      plot.subtitle = ggplot2::element_text(size = 10),
      panel.grid.minor = ggplot2::element_blank(),
      axis.title.y.right = ggplot2::element_text(margin = ggplot2::margin(l = 10)),
      legend.position = "top"
    )
}

plot_species_gg_onepanel("GRSP", target = 0.9)


#Estimate uncertainty around threshold ---------------------------------------------------------------------------------------------------------------------------

#Function for bootstrapping - NON-BI
bootstrap_threshold_nonbi <- function(species_name, n_boot = 100) {
  
  # Get species data
  sp_val <- val_bi_eval %>% filter(species == species_name)
  sp_labels <- lab_bins %>% filter(species == species_name)
  
  boot_results <- map_dfr(1:n_boot, function(b) {
    
    # Resample validation data
    boot_val <- sp_val %>% 
      slice_sample(n = nrow(sp_val), replace = TRUE)
    
    # Refit model
    boot_model <- tryCatch(
      glm(tp ~ score, family = binomial(), data = boot_val),
      error = function(e) NULL
    )
    
    if (is.null(boot_model)) {
      return(tibble(boot_iter = b, threshold = NA_real_, 
                    precision = NA_real_, prop_retained = NA_real_))
    }
    
    # Apply to FULL label bins (not resampled) - this is the population we're deploying to
    out <- compute_tseng_threshold(
      sp_labels %>% select(bin_lower, bin_mid, N_i),
      boot_model
    )
    
    tibble(
      boot_iter = b,
      threshold = out$summary$threshold,
      precision = out$summary$precision,
      prop_retained = out$summary$prop_retained
    )
  })
  
  # Summarize bootstrap distribution
  tibble(
    species = species_name,
    threshold_median = median(boot_results$threshold, na.rm = TRUE),
    threshold_ci_lower = quantile(boot_results$threshold, 0.025, na.rm = TRUE),
    threshold_ci_upper = quantile(boot_results$threshold, 0.975, na.rm = TRUE),
    precision_median = median(boot_results$precision, na.rm = TRUE),
    precision_ci_lower = quantile(boot_results$precision, 0.025, na.rm = TRUE),
    precision_ci_upper = quantile(boot_results$precision, 0.975, na.rm = TRUE),
    n_boot_success = sum(!is.na(boot_results$threshold))
  )
}

#run
bootstrap_nonbi <- map_dfr(unique(val_bi_eval$species), ~ bootstrap_threshold_nonbi(.x, n_boot = 100))
