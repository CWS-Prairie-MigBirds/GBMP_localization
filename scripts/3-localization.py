# -*- coding: utf-8 -*-
"""
by Erica Alex
2026-03-23

This script provides code for performing acoustic localization using the Opensoundscape library. Code is also provided for visualizing results, conducting error rejection, and preparing output results files.

Lapp, Sam; Rhinehart, Tessa; Freeland-Haynes, Louis; Khilnani, Jatin; Syunkova, Alexandra; Kitzes, Justin. “OpenSoundscape: An Open-Source Bioacoustics Analysis Package for Python.” Methods in Ecology and Evolution 2023. https://doi.org/10.1111/2041-210X.14196.

"""
import opensoundscape
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

#import csv with ARU coordinates ## UTM ZONE 13N

aru_coords = pd.read_csv("REDACTED", index_col=0)

#initialize a recorder array
from opensoundscape.localization import SynchronizedRecorderArray

array = SynchronizedRecorderArray(aru_coords)

#load csv of HawkEars detections
detections = pd.read_csv("REDACTED")

#Add the start timestamp of the recordings (should be the same for all receivers if the recordings are synchronized and trimmed)
import pytz
from datetime import datetime, timedelta

local_timestamp = datetime(2025, 6, 18, 6, 18, 59)
local_timezone = pytz.timezone("America/Regina")
detections["start_timestamp"] = [
    local_timezone.localize(local_timestamp) + timedelta(seconds=s)
    for s in detections["start_time"]
]

#set four column multi-index expected for localization
detections = detections.set_index(['file', 'start_time', 'end_time', 'start_timestamp'])

#set parameters for localization
min_n_receivers = 4 #min number of recievers with detection
max_receiver_dist = 80 #max distance between reference ARU and recievers used for localization

#localize detections
position_estimates5 = array.localize_detections(
    detections,
    min_n_receivers=min_n_receivers,
    max_receiver_dist=max_receiver_dist
)


# Data Exploration ==========================================================================================================

##Optionally filter for one species
#LCSP = [
#    e
#    for e in position_estimates5
#    if e.class_name == "LCSP"]

#Plot single position estimate and view details
example = position_estimates5[1]
print(f"The start time of the detection: {example.start_timestamp}")
print(f"This is a detection of the class/species: {example.class_name}")
print(f"The duration of the time-window in which the sound was detected: {example.duration}")
print(f"The estimated location of the sound: {example.location_estimate}")
print(f"The receivers on which our species was detected: \n{example.receiver_files}")
print(f"The estimated time-delays of arrival: \n{example.tdoas}")
print(f"The normalized Cross-Correlation scores: \n{example.cc_maxs}")

plt.plot(aru_coords["x"], aru_coords["y"], "^", label="ARU")
plt.scatter(
    x=example.location_estimate[0],
    y=example.location_estimate[1],
    color="red",
    label=f"{example.class_name}",
)
plt.legend(bbox_to_anchor=(1.02, 1), loc="upper left")
plt.show()

#Assess the quality of the localization based on how well the spectrogram lines up
from opensoundscape import Spectrogram

audio_segments = example.load_aligned_audio_segments()
specs = [Spectrogram.from_audio(a).bandpass(3000, 10000) for a in audio_segments]
plt.pcolormesh(np.vstack([s.spectrogram for s in specs]), cmap="Greys")

example.residual_rms #check residual

#See all the position estimates for this localized event
#each estimate is generated using a different recorder as the reference unit

sub = [
    e
    for e in position_estimates5
    if e.class_name == example.class_name
    and e.start_timestamp == example.start_timestamp
]

# get the x-coordinates of the estimated locations
x_coords = [e.location_estimate[0] for e in sub]
# get the y-coordinates of the estimated locations
y_coords = [e.location_estimate[1] for e in sub]
# get the rms of residuals per event
rms = [e.residual_rms for e in sub]
# plot the estimated locations, colored by the residuals
plt.scatter(
    x_coords,
    y_coords,
    c=rms,
    label="CCLO",
    alpha=0.4,
    edgecolors="black",
    cmap="jet",
)
cbar = plt.colorbar()
cbar.set_label("residual rms (meters)")
# plot the ARU locations
plt.plot(aru_coords["x"], aru_coords["y"], "^", label="ARU")
##Optionally plot known positions to compare
#known = pd.read_csv("REDACTED")
#plt.plot(168480.170614025,5448872.96142328, "X", markersize=5, label="Real position",c="red")
#plt.xlim(193700, 194050)
#plt.ylim(5.6216e6, 5.6220e6)
# make the legend appear outside of the plot
plt.legend(bbox_to_anchor=(1.2, 1), loc="upper left")

#Filter to only most precise localized positions based on residual rms
low_rms = [
    e for e in sub if e.residual_rms < 5
]  # get only the events with low residual rms

# plot the ARU locations
plt.plot(aru_coords["x"], aru_coords["y"], "^", label="ARU")
# plot the ARU locations
#plt.plot(-86.440, 58.861, "X", markersize=10, label="Real position")
# plot the estimated locations
plt.scatter(
    [e.location_estimate[0] for e in low_rms],
    [e.location_estimate[1] for e in low_rms],
    edgecolor="black",
    label="GRSP",
)
plt.legend(bbox_to_anchor=(1.02, 1), loc="upper left")

for e in low_rms:
    print(e.receiver_files)

#Plot all position estimates =============================================================================================================
x_coords = [e.location_estimate[0] for e in position_estimates5]
# get the y-coordinates of the estimated locations
y_coords = [e.location_estimate[1] for e in position_estimates5]
# get the rms of residuals per event
rms = [e.residual_rms for e in position_estimates5]
# plot the estimated locations, colored by the residuals
plt.scatter(
    x_coords,
    y_coords,
    c=rms,
    label="GRSP",
    alpha=0.4,
    edgecolors="black",
    cmap="jet",
)
cbar = plt.colorbar()
cbar.set_label("residual rms (meters)")
# plot the ARU locations
plt.plot(aru_coords["x"], aru_coords["y"], "^", label="ARU")
#plt.xlim(193745, 194015)
#plt.ylim(5621675,5621970)
# make the legend appear outside of the plot
plt.legend(bbox_to_anchor=(1.2, 1), loc="upper left")

type(position_estimates5)
dir(position_estimates5)
help(position_estimates5)
vars(position_estimates5)

#filter extreme outliers ==============================================================================================================
import numpy as np
from scipy.spatial import ConvexHull, distance

DIST_CUTOFF = 20  # buffer outside the convex hull of the array
RMS_CUTOFF = 20   # precision, NOTE: for playback and focal follows do not apply rms threshold

# Extract ARU coordinates
array_xy = aru_coords[["x", "y"]].values

# Compute convex hull
hull = ConvexHull(array_xy)
hull_vertices = array_xy[hull.vertices]

# Create a Path for quick inside/outside check
from matplotlib.path import Path
hull_path = Path(hull_vertices)

# Location estimates
loc_xy = np.array([e.location_estimate for e in position_estimates5])
rms = np.array([e.residual_rms for e in position_estimates5])

# Check if points are strictly inside hull
inside_hull = hull_path.contains_points(loc_xy)

# Compute distance from each point to the convex hull edges
def point_to_hull_distance(point, hull_vertices):
    """Return the minimum distance from point to edges of the hull."""
    min_dist = np.inf
    num_vertices = len(hull_vertices)
    for i in range(num_vertices):
        a = hull_vertices[i]
        b = hull_vertices[(i + 1) % num_vertices]
        # Distance from point to line segment ab
        seg_vec = b - a
        pt_vec = point - a
        t = np.clip(np.dot(pt_vec, seg_vec) / np.dot(seg_vec, seg_vec), 0, 1)
        proj = a + t * seg_vec
        dist = np.linalg.norm(point - proj)
        if dist < min_dist:
            min_dist = dist
    return min_dist

# Vectorized distance calculation
dist_to_hull = np.array([point_to_hull_distance(p, hull_vertices) for p in loc_xy])

# Keep points inside hull or within buffer distance
within_buffer = inside_hull | (dist_to_hull <= DIST_CUTOFF)

# Apply RMS filter
keep_rms = rms <= RMS_CUTOFF
keep = within_buffer & keep_rms

# Filtered positions
filtered_positions = [e for e, k in zip(position_estimates5, keep) if k]


#Check ===========================================================================================================
x_coords = [e.location_estimate[0] for e in filtered_positions]
# get the y-coordinates of the estimated locations
y_coords = [e.location_estimate[1] for e in filtered_positions]
# get the rms of residuals per event
rms = [e.residual_rms for e in filtered_positions]
#class_name
sp = [e.class_name for e in filtered_positions]
# Create mapping from class name → integer
unique_classes = list(set(sp))
class_to_int = {cls: i for i, cls in enumerate(unique_classes)}
sp_numeric = [class_to_int[s] for s in sp]

# plot the estimated locations, colored by the residuals
plt.scatter(
    x_coords,
    y_coords,
    c= rms,
    label="GRSP",
    alpha=0.4,
    edgecolors="black",
    cmap="jet",
)
cbar = plt.colorbar()
cbar.set_label("residual rms (meters)")
# plot the ARU locations
plt.plot(aru_coords["x"], aru_coords["y"], "^", label="ARU", c="yellow", mec="k", markersize=7)
#edit legend
plt.legend(bbox_to_anchor=(1.2, 1), loc="upper left")

#Playback =============================================================================================================================

#For playback recordings we will filter to only positions in which the reciever closest to the known position participated, and specify the start time of the playback. This ensures we are targeting the correct song event.  

final_positions = [
    e for e in filtered_positions
    if (
        (e.class_name == "WEME"
         and 39 in e.receiver_start_time_offsets
         and any("C2" in f for f in e.receiver_files))
        
        or
        
        (e.class_name == "VESP"
         and any(o in [54, 57] for o in e.receiver_start_time_offsets)
         and any("C2" in f for f in e.receiver_files))
        
        or
        
        (e.class_name == "CCSP"
         and any(o in [57, 60] for o in e.receiver_start_time_offsets)
         and any("C2" in f for f in e.receiver_files))
        
        or
        
        (e.class_name == "BRBL"
         and any(o in [63, 66] for o in e.receiver_start_time_offsets)
         and any("C2" in f for f in e.receiver_files))
        
        or
        
        (e.class_name == "LCSP"
         and 45 in e.receiver_start_time_offsets
         and any("C2" in f for f in e.receiver_files))
        
        or
        
        (e.class_name not in ["WEME", "VESP", "CCSP", "BRBL", "LCSP"])
    )
]


#For each species and start_timestamp, keep only the most precise position (lowest rms)
best_positions_ls = {}

for e in final_positions:
    key = (e.class_name, e.start_timestamp)
    
    # If key not seen OR this one has lower RMS → replace
    if key not in best_positions_ls or e.residual_rms < best_positions_ls[key].residual_rms:
        best_positions_ls[key] = e

# Final filtered list
best_positions = list(best_positions_ls.values())

#check plot
x_coords = [e.location_estimate[0] for e in best_positions]
# get the y-coordinates of the estimated locations
y_coords = [e.location_estimate[1] for e in best_positions]
# get the rms of residuals per event
rms = [e.residual_rms for e in best_positions]
#class_name
sp = [e.class_name for e in best_positions]
# Create mapping from class name → integer
unique_classes = list(set(sp))
class_to_int = {cls: i for i, cls in enumerate(unique_classes)}
sp_numeric = [class_to_int[s] for s in sp]

# plot the estimated locations, colored by the residuals
plt.scatter(
    x_coords,
    y_coords,
    c= sp_numeric,
    label="GRSP",
    alpha=0.4,
    edgecolors="black",
    cmap="jet",
)
cbar = plt.colorbar()
cbar.set_label("residual rms (meters)")
# plot the ARU locations
plt.plot(aru_coords["x"], aru_coords["y"], "^", label="ARU", c="yellow", mec="k", markersize=7)
#edit legend
plt.legend(bbox_to_anchor=(1.2, 1), loc="upper left")

# Follow ==============================================================================================================================

#For focal follow, filter to estimates in which the receiver closest to the known location participated

receiver_conditions = ["C2"]

final_positions = []

for e in filtered_positions:
    # check if any receiver condition is in any receiver file
    for f in e.receiver_files:
        if any(cond in f for cond in receiver_conditions):
            final_positions.append(e)
            break  # stop checking this e, already matched

#check plot
x_coords = [e.location_estimate[0] for e in final_positions]
# get the y-coordinates of the estimated locations
y_coords = [e.location_estimate[1] for e in final_positions]
# get the rms of residuals per event
rms = [e.residual_rms for e in final_positions]
#class_name
sp = [e.class_name for e in final_positions]
# Create mapping from class name → integer
unique_classes = list(set(sp))
class_to_int = {cls: i for i, cls in enumerate(unique_classes)}
sp_numeric = [class_to_int[s] for s in sp]

# plot the estimated locations, colored by the residuals
plt.scatter(
    x_coords,
    y_coords,
    c= rms,
    label="GRSP",
    alpha=0.4,
    edgecolors="black",
    cmap="jet",
)
cbar = plt.colorbar()
cbar.set_label("residual rms (meters)")
# plot the ARU locations
plt.plot(aru_coords["x"], aru_coords["y"], "^", label="ARU", c="yellow", mec="k", markersize=7)
#edit legend
plt.legend(bbox_to_anchor=(1.2, 1), loc="upper left")


#To assess localization accuracy, we will match the localized positions to observations made during the focal follow based on closest matching timestamps
from datetime import datetime
import pandas as pd

# load csv of observations and modify date/time to match localized positions
follow = pd.read_csv("REDACTED")
follow['start_timestamp'] = pd.to_datetime(follow['date'] + ' ' + follow['time'])
csv_timestamps = follow['start_timestamp'].tolist()

# make sure position_estimates have datetime objects
for e in final_positions:
    if isinstance(e.start_timestamp, str):
        e.start_timestamp = datetime.fromisoformat(e.start_timestamp)

# align timezones
position_tz = final_positions[0].start_timestamp.tzinfo
csv_timestamps_tz = [
    ts.tz_localize(position_tz) if ts.tzinfo is None else ts.astimezone(position_tz)
    for ts in csv_timestamps
]

# find nearest matching timestamps, keep the position_estimate with the lowest residual rms
best_positions = []
for ts in csv_timestamps_tz:
    # Find positions with minimum absolute difference in start_timestamp
    min_diff = min(abs(e.start_timestamp - ts) for e in final_positions)
    closest_candidates = [e for e in final_positions if abs(e.start_timestamp - ts) == min_diff]
    
    # Pick the one with the lowest residual_rms among candidates
    best_candidate = min(closest_candidates, key=lambda e: e.residual_rms)
    best_positions.append(best_candidate)


print(f"Number of CSV timestamps: {len(csv_timestamps_tz)}")
print(f"Number of final matched positions: {len(best_positions)}")  # should match CSV

#Note: Some song events may not have been localized successfully, so two distinct observations may match to the same positon_estimate. Make sure to remove duplicates as necessary.

#Check plot
x_coords = [e.location_estimate[0] for e in best_positions]
# get the y-coordinates of the estimated locations
y_coords = [e.location_estimate[1] for e in best_positions]
# get the rms of residuals per event
rms = [e.residual_rms for e in best_positions]
#class_name
sp = [e.class_name for e in best_positions]
# Create mapping from class name → integer
unique_classes = list(set(sp))
class_to_int = {cls: i for i, cls in enumerate(unique_classes)}
sp_numeric = [class_to_int[s] for s in sp]

# plot the estimated locations, colored by the residuals
plt.scatter(
    x_coords,
    y_coords,
    c= rms,
    label="GRSP",
    alpha=0.4,
    edgecolors="black",
    cmap="jet",
)
cbar = plt.colorbar()
cbar.set_label("residual rms (meters)")
# plot the ARU locations
plt.plot(aru_coords["x"], aru_coords["y"], "^", label="ARU", c="yellow", mec="k", markersize=7)
#known positions
plt.plot(193818.1, 5621742, marker="X", markersize=8, c="red", label = "Known Position")
#edit legend
plt.legend(bbox_to_anchor=(1.2, 1), loc="upper left")

# Mini spec review ====================================================================================================================

#We will conduct minimum spectrogram review for all general localizations (e.g. for full 30min recs, not focal follows or playback) as a quality control measure. This process collects audio from all of the recievers that participated in #a given localization event and converts them to spectrograms. For each pixel position, the lowest value across recievers is selected and compiled into a single "minimum spectrogram" (minspec) for that event. The minspecs are then #analyzed by HawkEars and filtered based on a score threshold. Events where the species detection was incorrect, or audio cross-correlation failed (poor alignment), the score HawkEars score will be low, and subsequently filtered out.

import numpy as np
import librosa
from pathlib import Path as PLPath  # avoid conflict with matplotlib.path.Path
from joblib import Parallel, delayed
from opensoundscape import Audio, Spectrogram

#Parameters
temp_audio_dir = "REDACTED" #path to temporarily save minspecs
minspec_discard_distance = 45  #max distance between recievers
n_jobs = 12  # parallel workers

#Helper Functions
def spec_to_audio(spec, sr):
    """Convert a spectrogram back to Audio using Griffin-Lim."""
    y_inv = librosa.griffinlim(spec.spectrogram, hop_length=256, win_length=512)
    return Audio(y_inv, sr)

def distances_to_receivers(position, dims=2):
    """Compute distances from the position estimate to each receiver."""
    return [
        np.linalg.norm(position.location_estimate[:dims] - r[:dims])
        for r in position.receiver_locations
    ]

def min_spec_to_audio(position, discard_over_distance=50):
    """Generate a minimum spectrogram audio clip from aligned segments."""
    clips = position.load_aligned_audio_segments()
    distances = distances_to_receivers(position)

    # Keep only nearby receivers
    clips = [c for i, c in enumerate(clips) if distances[i] < discard_over_distance]

    # If no valid clips → fallback
    if len(clips) == 0:
        raise ValueError("No receivers within distance")

    # Convert clips to spectrograms
    specs = [Spectrogram.from_audio(c, dB_scale=False) for c in clips]

    # Minimum spectrogram across clips
    minspec = specs[0]._spawn(
        spectrogram=np.min(np.array([s.spectrogram for s in specs]), axis=0)
    )

    # Normalize to loudest clip
    max_val = np.max([c.samples.max() for c in clips])

    return (
        spec_to_audio(minspec, clips[0].sample_rate)
        .normalize(max_val)
        .extend_to(clips[0].duration)
    )

#Create temp directory
clip_dir = PLPath(temp_audio_dir)
clip_dir.mkdir(parents=True, exist_ok=True)

print(f"Generating minspec clips for {len(filtered_positions)} positions...")

#Process each position in parallel
def process_position(position, idx):
    """Generate min-spec audio for a position, or fallback to silence."""
    try:
        out_path = clip_dir / f"{idx}.wav"
        min_spec_to_audio(position, discard_over_distance=minspec_discard_distance).save(out_path)
        return 0
    except Exception as e:
        # fallback to 3 seconds silence
        Audio.silence(sample_rate=44100, duration=3).save(clip_dir / f"{idx}.wav")
        return 1

# Run in parallel
results = Parallel(n_jobs=n_jobs)(
    delayed(process_position)(p, i) for i, p in enumerate(filtered_positions)
)

print(f"Failures: {sum(results)} / {len(results)}")

#Run HawkEars over the minspecs, then come back =======================================================================================================================================================================

#Pull in HawkEars detections
hawkears_results = pd.read_csv("REDACTED")

hawkears_results["clip_id"] = (
    hawkears_results["filename"]
    .str.replace(".wav", "", regex=False)
    .astype(int)
)

scores = hawkears_results[['clip_id', 'score']]

# df_scores: columns ['clip_id', 'score']
score_dict = dict(zip(scores["clip_id"], scores["score"]))

#The ID of the minspec corresponds to it's index position in the localization data. Use this to join scores to position estimates
for idx, pos in enumerate(filtered_positions):
    pos.score = score_dict.get(idx, None)  # None if missing


import numpy as np
from scipy.spatial.distance import cdist
from collections import defaultdict

#set thrshold for HawkEars score
SCORE_THRESHOLD = 0.4
filtered_by_score = [
    pos for pos in filtered_positions
    if pos.score is not None and pos.score >= SCORE_THRESHOLD
]

#We will also aggregate localized positions that are assumed to be the same individual. For localized events that are within 10m of each other, and have the same class name and start timestamp, we will keep only the position with the lowest residual rms. This will reduce the number of duplicate localizations for the same bird, while ensuring that spatially explicit song events are preserved.

DEDUP_DISTANCE = 10  # meters

# Group positions by start_timestamp and class_name
timestamp_groups = defaultdict(list)
for pos in filtered_by_score:
    key = (pos.start_timestamp, pos.class_name)
    timestamp_groups[key].append(pos)

deduped_positions = []

for ts, positions in timestamp_groups.items():
    if len(positions) == 1:
        deduped_positions.append(positions[0])
        continue
    
    # Create array of positions' coordinates
    locs = np.array([p.location_estimate for p in positions])
    
    # Flags to track which positions to keep
    keep_flags = np.ones(len(positions), dtype=bool)
    
    # Compute pairwise distances
    dist_matrix = cdist(locs, locs)
    
    for i in range(len(positions)):
        if not keep_flags[i]:
            continue
        for j in range(i + 1, len(positions)):
            if not keep_flags[j]:
                continue
            if dist_matrix[i, j] <= DEDUP_DISTANCE:
                # Keep the one with the lower residual_rms
                if positions[i].residual_rms <= positions[j].residual_rms:
                    keep_flags[j] = False
                else:
                    keep_flags[i] = False
    
    # Append surviving positions
    for i, keep in enumerate(keep_flags):
        if keep:
            deduped_positions.append(positions[i])

# Replace filtered_positions with final filtered/deduplicated list
final_positions = deduped_positions

#Check plot
print(f"Number of positions after filtering: {len(final_positions)}")

#plot
x_coords = [e.location_estimate[0] for e in final_positions]
# get the y-coordinates of the estimated locations
y_coords = [e.location_estimate[1] for e in final_positions]
# get the rms of residuals per event
rms = [e.residual_rms for e in final_positions]
# plot the estimated locations, colored by the residuals
plt.scatter(
    x_coords,
    y_coords,
    c=rms,
    label="CCLO",
    alpha=0.4,
    edgecolors="black",
    cmap="jet",
)
cbar = plt.colorbar()
cbar.set_label("residual rms (meters)")
# plot the ARU locations
plt.plot(aru_coords["x"], aru_coords["y"], "^", label="ARU")
# make the legend appear outside of the plot
plt.legend(bbox_to_anchor=(1.2, 1), loc="upper left")

# Convert PositionEstimate objects to a DataFrame and save audio clips =======================================================================================

import pandas as pd
import os
from opensoundscape import Audio

#Parameters - save 10s audio clips centred around the 3s song event
CLIP_DURATION = 10.0
HALF_WINDOW = CLIP_DURATION / 2
DETECTION_WINDOW = 3.0
OUTPUT_DIR = "REDACTED" #where to save audio clips
AUDIO_FOLDER_NAME = "audio"
os.makedirs(OUTPUT_DIR, exist_ok=True)

#build dataframe from positions object
rows = []
for e in best_positions:  # final_positions should be a list of PositionEstimate-like objects
    rows.append({
        "class_name": e.class_name,
        "start_timestamp": e.start_timestamp,
        "duration": e.duration,
        "position": list(e.location_estimate),
        "receiver_files": list(e.receiver_files),
        "receiver_start_time_offset": list(e.receiver_start_time_offsets),
        "cc_maxs": list(e.cc_maxs),
        "tdoas": list(e.tdoas),
        "residual_rms": e.residual_rms
        #"hawkears_score": e.score
    })

filtered_df = pd.DataFrame(rows)

##Create a unique event ID to match localized positions to audio clips etc

#for general 30min localizations, append classname and datetime info
#base_id = (
#    filtered_df["class_name"].astype(str) + "_",
#    filtered_df["start_timestamp"].dt.strftime("%Y%m%dT%H%M%S")
#)
#filtered_df["event_id"] = base_id + "_" + filtered_df.groupby(base_id).cumcount().add(1).astype(str)

#for playback and follows, append observation type
base_id = (
    filtered_df["class_name"].astype(str)
    + "_follow_"
    + filtered_df["start_timestamp"].dt.strftime("%Y%m%dT%H%M%S")
)

filtered_df["event_id"] = (
    base_id
    + "_"
    + filtered_df.groupby(base_id).cumcount().add(1).astype(str)
)

# Reorder columns
filtered_df = filtered_df[
    [
        "event_id", "class_name", "start_timestamp", "duration",
        "position", "receiver_files",
        "receiver_start_time_offset", "cc_maxs", "tdoas", "residual_rms" #"hawkears_score"
    ]
]

#Clip audio windows for each localization event, save to folders that correspond to ARU IDs. Map relative path and append to dataframe
clip_records = []
receiver_to_clip = {}

# Lists to store updated columns for filtered_df
updated_receiver_files = []
updated_receiver_offsets = []

for idx, row in filtered_df.iterrows():
    event_id = row["event_id"]
    start_timestamp = row["start_timestamp"]
    receiver_files = row["receiver_files"]
    receiver_offsets = row["receiver_start_time_offset"]

    saved_clips_for_row = []
    adjusted_offsets_for_row = []

    for rec_path, offset in zip(receiver_files, receiver_offsets):
        rec_path = os.path.normpath(rec_path)
        if not os.path.exists(rec_path):
            print(f"Skipping {rec_path}: file does not exist.")
            continue

        clip_start = float(offset) - HALF_WINDOW

        try:
            audio = Audio.from_file(rec_path, offset=clip_start, duration=CLIP_DURATION)

            rec_name = os.path.splitext(os.path.basename(rec_path))[0]
            aru_id = rec_name.split("_")[0]

            # ARU subfolder
            aru_folder = os.path.join(OUTPUT_DIR, aru_id)
            os.makedirs(aru_folder, exist_ok=True)

            clip_filename = f"{event_id}_{rec_name}.flac"
            out_path = os.path.join(aru_folder, clip_filename)
            audio.save(out_path)

            # Save just the filename
            saved_clips_for_row.append(clip_filename)
            receiver_to_clip[rec_path] = clip_filename

            # Offset relative to clip start
            relative_offset = HALF_WINDOW - (DETECTION_WINDOW / 2)
            adjusted_offsets_for_row.append(relative_offset)

            # Clip index CSV
            rel_path_from_audio = os.path.join(AUDIO_FOLDER_NAME, os.path.relpath(out_path, OUTPUT_DIR))
            clip_records.append({
                "filename": clip_filename,
                "relative_path": rel_path_from_audio,
                "recorder_id": aru_id,
                "start_timestamp": start_timestamp
            })

        except Exception as e:
            print(f"Failed to clip {rec_path} for {event_id}: {e}")

    # Append updated info
    updated_receiver_files.append(saved_clips_for_row)
    updated_receiver_offsets.append(adjusted_offsets_for_row)

#updated dataframe with new columns - reciever files and start_time_offsets relative to audio clips rather than raw audio files
filtered_df["receiver_files"] = updated_receiver_files
filtered_df["receiver_start_time_offset"] = updated_receiver_offsets

#save audio clip index
clips_df = pd.DataFrame(clip_records)
clip_index_csv = os.path.join(OUTPUT_DIR, "REDACTED")
clips_df.to_csv(clip_index_csv, index=False)
print(f"Saved clip index CSV: {clip_index_csv}")

#save csv of final, filtered localized positions
output_csv = "REDACTED"
filtered_df.to_csv(output_csv, index=False)
print(f"Saved filtered positions CSV: {output_csv}")

## Optionally plot the spectrograms to visualize alignment ================================================================================================================================================

import os
import matplotlib.pyplot as plt
from opensoundscape import Spectrogram
from tqdm import tqdm
import numpy as np

SAVE_FIGS = False  #True to save, False to only view
OUTPUT_DIR = "REDACTED"
os.makedirs(OUTPUT_DIR, exist_ok=True)

def get_recorder_name(file_path):
    """Extract recorder prefix from file path"""
    base = os.path.basename(file_path)
    name = base.split("_")[0]
    return name

# Loop over filtered events from CSV
for idx, row in tqdm(filtered_df.iterrows(), total=len(filtered_df)):

    event_id = row["event_id"]

    # Select PositionEstimate matching class, timestamp, AND RMS
    pe = next(
        (e for e in final_positions
         if e.class_name == row["class_name"]
         and e.start_timestamp == row["start_timestamp"]),
        None
    )
    if pe is None:
        print(f"Warning: PositionEstimate not found for {event_id}")
        continue

    # Load aligned audio segments
    try:
        audio_segments = pe.load_aligned_audio_segments()
    except Exception as e:
        print(f"Failed to load audio segments for {event_id}: {e}")
        continue

    # Compute recorder names
    receiver_names = [get_recorder_name(f) for f in pe.receiver_files]

    # Compute spectrograms
    specs = []
    names_valid = []
    for seg, name in zip(audio_segments, receiver_names):
        try:
            spec = Spectrogram.from_audio(seg).bandpass(6000, 12000)
            if spec.spectrogram.shape[1] < 1:
                continue  # skip empty segments
            specs.append(spec)
            names_valid.append(name)
        except Exception as e:
            print(f"Skipping receiver {name} for event {event_id}: {e}")
            continue

    if len(specs) == 0:
        print(f"No valid receivers for {event_id}, skipping plot")
        continue

    # Stack spectrograms vertically
    stacked_spec = np.vstack([s.spectrogram for s in specs])

    # Plot
    plt.figure(figsize=(12, 6))
    plt.pcolormesh(stacked_spec, cmap="Greys")
    plt.title(f"Aligned Spectrograms: {event_id}| Residual RMS: {pe.residual_rms:.2f} m| HawkEars Score: {pe.score:.2f}")
    plt.xlabel("Time bins")
    plt.ylabel("Receivers")

    # Label y-axis with recorder names
    n_receivers = len(names_valid)
    plt.yticks(
        ticks=np.linspace(0, stacked_spec.shape[0], n_receivers, endpoint=False) 
              + stacked_spec.shape[0]/(2*n_receivers),
        labels=names_valid
    )

    plt.tight_layout()

    # Save or show
    if SAVE_FIGS:
        save_path = os.path.join(OUTPUT_DIR, f"{event_id}_spectrogram.png")
        plt.savefig(save_path, dpi=150)
        plt.close()
    else:
        plt.show()
