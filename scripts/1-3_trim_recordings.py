# -*- coding: utf-8 -*-
"""
by Tessa Rhinehart
2024-11-11
https://github.com/rhine3/frontierlabsutils

Modified by Erica Alex
2026-02-10

This script trims resampled recordings (1-2_sync_recordings.py) and trims them so that they contain audio from the same time window. This is accomplished by finding the latest start time and earliest end time across recordings, and trimming all recordings to that time.

Important: modify and run frontierlabsutils.py before running this script

"""

import pandas as pd
from pathlib import Path
from opensoundscape.audio import Audio

from frontierlabsutils import get_recording_path, get_recorder_list, get_all_times
from frontierlabsutils import get_latest_start_second, get_earliest_end_second
from frontierlabsutils import extract_start_end, get_audio_from_time

data_dir =  "REDACTED" # Directory containing resampled files
out_dir =  "REDACTED" # Directory to save trimmed recordings

# Make the trimmed recording directory
trimmed_recording_dir = Path(out_dir)
trimmed_recording_dir.mkdir(exist_ok=True)

#Specify date and times to find desired recordings
recorders = get_recorder_list()
date = "20250620"
earliest_start_times = ["0617"]

#look for recordings within 3s of the specified times, can expand this value if drift was more severe
for start_time in earliest_start_times:
    try_times = get_all_times(start_time)[:3]
    recordings_this_time = []
    recorders_this_time = []
    for recorder in recorders:
        for time in try_times:
            recording = get_recording_path(
                recorder=recorder,
                date=date,
                hour_minute=time,
                data_dir=data_dir
            )
            if recording:
                recorders_this_time.append(recorder)
                recordings_this_time.append(Path(recording))
                break

    # Find the latest start time across recordings
    latest_start_second = get_latest_start_second(recordings_this_time)

    # Find the earliest end time across recordings
    earliest_end_second = get_earliest_end_second(recordings_this_time)
    
    print(f"Latest start time across recorders: {latest_start_second}")
    print(f"Earliest end time across recorders: {earliest_end_second}")

    # Trim all recordings from the earliest time to the end
for recorder, recording in zip(recorders_this_time, recordings_this_time):
    
    # Create folder: trimmed_recording_dir / date / earliest_start_time
    folder_name = f"{date}_{latest_start_second.strftime('%H%M')}"
    recorder_out_dir = trimmed_recording_dir / folder_name
    recorder_out_dir.mkdir(parents=True, exist_ok=True)  # make parent folders too
    
    # Add recorder name as prefix to file
    new_filename = f"{recorder}_{recording.name}"
    trimmed_audio_filename = recorder_out_dir / new_filename

    if trimmed_audio_filename.exists():
        print(f"Already trimmed: {trimmed_audio_filename.name}")
        continue

    print(f"Trimming {recording.name}...")

    audio = Audio.from_file(recording)
    
    original_start, original_end = extract_start_end(recording.name)
    
    clip_len = (earliest_end_second - latest_start_second).seconds
    
    trimmed_audio = get_audio_from_time(
        clip_start = latest_start_second,
        clip_length_s = clip_len,
        original_start = original_start,
        original_audio = audio
    )
    trimmed_audio.save(trimmed_audio_filename)
    print(f"Saved trimmed audio: {trimmed_audio_filename}")