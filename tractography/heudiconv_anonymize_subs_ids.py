#!/usr/bin/env python
import sys
import os
import pandas as pd

INPUT_DATA_FOLDER = os.path.join("data", "original_dicom")

subjs_names = next(os.walk(INPUT_DATA_FOLDER))[1]
num_subjs = len(subjs_names)
subjs_ids = list(range(num_subjs))
subjs_ids_str = [str(id) for id in subjs_ids]

subj_map = dict(zip(subjs_names, subjs_ids_str))
subj_map_pd = pd.DataFrame([subjs_names, subjs_ids_str]).T
subj_map_pd.columns = ["name", "id"]
subjs_map_filepath = os.path.join("data", "heudiconv_subjs_anonymized_map.csv")
subj_map_pd.to_csv(subjs_map_filepath, index=False)

sid = sys.argv[-1]
if sid in subj_map:
    print(subj_map[sid])