#Creating csv for Opensoundscape from HawkEars detections

#by Erica Alex
#2025-12-06

#This script prepares csv files in the appropriate format for conducting acoustic localization using the Opensoundscape python library
#https://opensoundscape.org/en/latest/tutorials/acoustic_localization.html

#The library expects a detections.csv with the following columns:
#file: path to audio file for each receiver
#start_time: start of time window (in seconds)
#end_time: end of time window (in seconds) ***Note: the length of time windows should be consistent, and should encompass the entire recording length
#start_timestamp: start time of the recording in ISO format (should be the same across all recorders if synchronized and trimmed). Can be created in the localization.py script
#species classes: one column for every species of interest. Set of binary detections for each species in each time window.

#and an aru_coords.csv with the following columns:
#file: path to audio file for each reciever
#x: x position of the receiver (UTM)
#y: y position of the reciever (UTM)

library(tidyverse)
library(stringr)
library(sf)

#Transform HawkEars detections to binary csv for localization --------------------------------------------------------------------------------------

labels <- read.csv("REDACTED")  #csv with HawkEars detections

#str(labels)

#optionally filter for particular species
#playback <- c("NSWO", "DEJU", "AMRO", "WEME", "MOWA", "AMCR", "LCSP", "YEWA", "BADO", "ALFL", "VESP", "CCSP",
#              "BOCH", "BRBL", "MAWR")

#labels = labels %>%
#   filter(class_code %in% playback)

#create tbl for all time windows
time_grid <- tibble(
  start_time = seq(0, 1629, by = 3),
  end_time   = start_time + 3
)

#get all filenames
files <- labels %>%
  distinct(filename)

#create full filename × time grid
full_grid <- files %>%
  crossing(time_grid)

#get presence/absence from HE detections
presence <- labels %>%
  distinct(filename, start_time, end_time, class_code) %>%
  mutate(present = 1)

#join and pivot to wide 
final_df <- full_grid %>%
  left_join(presence,
            by = c("filename", "start_time", "end_time")) %>%
  mutate(present = if_else(is.na(present), 0L, present)) %>%
  pivot_wider(
    names_from  = class_code,
    values_from = present,
    values_fill = 0
  ) %>%
  arrange(filename, start_time)

#join filepaths from basename
files <- list.files("REDACTED", full.names = TRUE)

final_df$file <- files[match(final_df$filename, basename(files))]

final_df = final_df %>%
  select(file, start_time, end_time, CCLO)

#save
write.csv(final_df, "REDACTED", row.names=FALSE)

#Create ARU coords ----------------------------------------------------------------------------------------------------------------------------------

coords <- read.csv("REDACTED") #read in list of receivers and their positions (lat/lon from RTK)

#Convert lat/lon to UTM

coords_sf <- st_as_sf(
  coords,
  coords = c("longitude", "latitude"),
  crs = 4326)

#UTM Zone 13N
coords_utm13 <- st_transform(coords_sf, crs = 32613)

#extract to columns
coords_utm13 <- cbind(
  coords_utm13,
  st_coordinates(coords_utm13)
)

#drop geometry
coords_utm13_df <- st_drop_geometry(coords_utm13)

#join filepaths based on aru_id 
coords <- coords %>%
  rowwise() %>%
  mutate(file = list(files[str_starts(basename(files), aru_id)])) %>%
  unnest(cols = c(file))

coords = coords %>%
  select(file, x, y)

#save file
write.csv(coords, "REDACTED", row.names = FALSE)
