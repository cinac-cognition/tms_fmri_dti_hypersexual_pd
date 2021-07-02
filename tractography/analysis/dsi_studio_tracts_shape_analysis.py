import subprocess
from pathlib import Path
from shutil import copyfile


# import numpy as np
import os
# import pandas as pd

print(" ================= PATH ===============")
print(os.getcwd())

# derivatives_data = os.path.join("..", "data", "derivatives")  # For local machine
derivatives_data = os.path.join("data", "data", "derivatives")  # For Docker
dwi_data = os.path.join(derivatives_data, "qsiprep")
dti_data = os.path.join(derivatives_data, "dti")
tracts_data = os.path.join(derivatives_data, "tracts")
# The "hcp1065_2mm.fib.gz" is obtained from the dsi studio installation folder
atlas_fib_path = os.path.join(tracts_data, "atlas", "hcp1065_2mm.fib.gz")

if not os.path.exists(tracts_data):
    os.mkdir(tracts_data)

num_subjs = 0 # 0 = ALL
subjs_names = sorted([f.name for f in Path(tracts_data).iterdir() if f.is_dir() and f.name != "atlas"])
if num_subjs != 0:
    subjs_names = subjs_names[0:num_subjs]

def get_data(subj, name, type):
    path_dir = ""
    if type == "dwi":
        path_dir = os.path.join(dwi_data, subj, "ses-0", "dwi")
    elif type == "dti":
        path_dir = os.path.join(dti_data, subj)
    elif type == "tracts":
        path_dir = os.path.join(tracts_data, subj)

    path = ""

    if name == "whole_tract":
        path = "{}_whole_trk.trk".format(subj)
    elif name == "whole_tract_mni":
        path = "{}_whole_trk_mni.trk".format(subj)
    elif name == "rec_bundles_dir":
        path = "rec_bundles"
    elif name == "rec_bundles":
        path = os.path.join("rec_bundles", "*.trk")
    elif name.startswith("tract__"):
        tract_name = name.split("tract__")[1]
        path = os.path.join("org_bundles", tract_name)
    elif name.startswith("atlas_tract__"):
        tract_name = name.split("tract__")[1]
        tract_filename = "{}.trk".format(tract_name)
        path = os.path.join("rec_bundles", tract_filename)

    path = os.path.join(path_dir, path)

    return path

def run_cmd(dipy_full_cmd):
    global p, line
    p = subprocess.Popen(dipy_full_cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)

    while True:
        line = p.stdout.readline()
        print(">>>", line)
        if not line: break


def bundles_shape_analysis():
    # http://dsi-studio.labsolver.org/Manual/command-line-for-dsi-studio#TOC-Tract-specific-analysis-voxel-based-analysis-connectivity-matrix-and-network-measures
    global i, subj, subj_path, output_dir, p, line

    for i, subj in enumerate(subjs_names):
        print("----- bundles_shape_analysis for subj {} ({} of {})".format(subj, i + 1, len(subjs_names)))
        subj_path = os.path.join(tracts_data, subj)

        subj_bundles_dir = get_data(subj, "rec_bundles_dir", type="tracts")
        tracts_names = sorted([f.name for f in Path(subj_bundles_dir).iterdir() if f.suffix == ".trk"])

        for j, tract_name in enumerate(tracts_names):
            print("Generating shape stats for tract {} ({} of {})".format(tract_name, j + 1, len(tracts_names)))
            tract_name = tract_name.split("__")[0] + "__labels__recognized_orig.trk"
            tract_path = get_data(subj, "tract__{}".format(tract_name), type="tracts")

            dsi_studio_cmd = "dsi_studio --action=ana --source={} --tract={} --export=stat".format(
                atlas_fib_path, tract_path
            ).split(" ")

            run_cmd(dsi_studio_cmd)

    return 0

if __name__ == "__main__":
    bundles_shape_analysis()