import json
import os

def set_all_intended_for():
    task = "task_alien_cookies"
    root_data_path = "data"
    data_path = os.path.join(root_data_path, "bids")
    data_files_placeholder = [
        data_path + os.sep + "sub-{}" + os.sep + "ses-0" + os.sep + "fmap" + \
                            os.sep + "sub-{}_ses-0_dir-AP_epi.json",
        data_path + os.sep + "sub-{}" + os.sep + "ses-0" + os.sep + "fmap" + \
        os.sep + "sub-{}_ses-0_dir-PA_epi.json"
        ]


    subjs_ids = [fpath.name.split("-")[1] for fpath in os.scandir(data_path) if fpath.is_dir()]

    for subj in subjs_ids:
        print("-------- Setting intended for in subj ", subj)

        for data_file in data_files_placeholder:
            fmap_json_path = data_file.format(subj, subj)
            with open(fmap_json_path, "r") as json_file:
                fmap_json = json.load(json_file)

            intended_for_vals = []
            dwi_run_path = "ses-0/dwi/sub-{}_ses-0_run-01_dwi.nii.gz".format(subj)
            intended_for_vals.append(dwi_run_path)

            fmap_json["IntendedFor"] = intended_for_vals
            with open(fmap_json_path, "w") as json_file:
                json.dump(fmap_json, json_file)


if __name__ == "__main__":
    print("Starting setting fieldmap intended for field in json sidecars...")
    set_all_intended_for()
    print("Done setting all intended for")