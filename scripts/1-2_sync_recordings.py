# -*- coding: utf-8 -*-
"""
by Tessa Rhinehart
2024-11-11
https://github.com/rhine3/frontierlabsutils

Modified by Erica Alex
2026-02-02

Important: modify and run the frontierlabsutils.py script first

This script provides code for synchronizing BAR-LT recordings for acoustic localization by identifying buffer overflows and dropped buffers, and inserting blank space (adding zeros).

"""
from opensoundscape.audio import Audio
from pathlib import Path
import frontierlabsutils
print(frontierlabsutils.__file__)
from time import time
import soundfile as sf

data_dir = "REDACTED" # Where your original recordings and loclog files are stored
out_dir  = "REDACTED"# Where to save your resampled recordings (created below)

# Define recorders A1 to G7, optionally add a prefix for array ID. Must match frontierlabsutils.py
recorders = [f"SNAS-{letter}{num}" for letter in "ABCDEFG" for num in range(1, 8)]

# Dates you want to synchronize recordings from, optionally specify a time (or synchronize all recs for that date)
dates = ["20250627"]
#times = ["0618"]

resampled_folder = Path(out_dir)

t0 = time()
num_resampled = 0

for recorder in recorders:
    print(f"\nProcessing recorder: {recorder}")

    # Make sure the output folder exists for this recorder
    recorder_folder = resampled_folder.joinpath(recorder)
    recorder_folder.mkdir(parents=True, exist_ok=True)

    for date in dates:
        print(f"\nResampling data from recorder {recorder}, date: {date}")

        # Get all recordings from this date and recorder
        rec_date_recording_paths = frontierlabsutils.get_recording_path(
            recorder=recorder, 
            date=date,
            hour_minute="any", 
            data_dir=data_dir
        )

        # Avoid dates with zero recordings found
        if type(rec_date_recording_paths) != list or len(rec_date_recording_paths) == 0:
            print(f"  No recordings found for {recorder} on {date}. Skipping.")
            continue
            
        print(f"{len(rec_date_recording_paths)} recordings found")

        # Get loclog contents
        loclog_values = frontierlabsutils.get_loclog_contents(
            data_dir=data_dir,
            recorder=recorder,
            date=date
        )
        
        if len(loclog_values) < 1:
            print("  No loclog values found for recorder. Continuing.")
            continue

        # Resample all recordings
        for recording_path in rec_date_recording_paths:
            print(f"Resampling {recording_path}")
            if sf.info(recording_path).duration > 3600: print("Skipping corrupt WAV"); continue

            resampled_filename = recorder_folder.joinpath(f"{Path(recording_path).stem}_resampled.wav")
            start, end = frontierlabsutils.extract_start_end(Path(recording_path).name)
            

            if resampled_filename.exists():
                print("  Already resampled. Continuing.")
                continue

            try:
                recording_time = recording_path.name.split('.')[0].split('T')[1][:4]
            except IndexError:
                print(f"  Recording {recording_path} not in correct format. Skipping.")
                continue

            write_times = frontierlabsutils.get_recording_write_times(
                loclog_values, date=date, start_time=recording_time)
            if len(write_times) == 0:
                print(f"  No write times logged for {recorder}. Skipping.")
                continue

            overflow_sample_indices, overflow_samples_to_insert = frontierlabsutils.get_buffer_insert_lengths(write_times)            
            audio = frontierlabsutils.insert_missing_buffers(
                recorder=recorder,
                date=date,
                time=recording_time,
                data_dir=data_dir,
                overflow_sample_indices=overflow_sample_indices,
                overflow_samples_to_insert=overflow_samples_to_insert
            )

            resampled = frontierlabsutils.interpolate_audio(start, end, audio)
            resampled.save(resampled_filename)
            num_resampled += 1

t1 = time()
print(f"\nIt took {(t1-t0)/60:.2f} minutes to resample {num_resampled} files.")
