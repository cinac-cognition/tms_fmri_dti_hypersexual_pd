#!/bin/bash

DATA_INPUT_FOLDER="data/original_dicom"
DATA_OUTPUT_FOLDER="data/bids"

\rm -R $DATA_OUTPUT_FOLDER
\rm -R .heudiconv

SUBJECTS_LIST=$(find $DATA_INPUT_FOLDER -maxdepth 1 -mindepth 1 -type d |  cut -f3 -d'/')
echo $SUBJECTS_LIST
echo "----- Make sure you DO NOT have the DICOM series data folders in SUBFOLDERS inside every subject folder -------"

heudiconv \
-d $DATA_INPUT_FOLDER/{subject}/*/*.dcm \
-s $SUBJECTS_LIST \
-ss 0 \
-f ./heudiconv_config.py \
-c dcm2niix \
-o $DATA_OUTPUT_FOLDER \
-b --overwrite \
--minmeta

# Consider including the following argument in casse you want to anonymize the data:
# --anon-cmd ./heudiconv_anonymize_subs_ids.py \

\rm -R $DATA_OUTPUT_FOLDER/.heudiconv

python heudiconv_set_fieldmap_intended_for.py

# If it doesn't work, try running:
# python fix_epi_different_sessions

# Without --minmeta it can fail when running fmriprep SDC with fieldmap

# python convert_task_events_csv.py
