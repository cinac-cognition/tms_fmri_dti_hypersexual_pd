## Project description
DTI for PD Hypersexual vs PD Non-Hypersexual patients.
Project for PhD UAM-HM CINAC.
*Only scripts are included in this repository. Ask the authors for the data if needed.

DTI pipeline:
- DWI preprocessing (denoise, correct Eddy currents, etc): qsiprep
- DTI estimation (with seeds, end conditions, forbidden pahts, etc): dipy NLLS
- DTI postprocessing: tensor normalization with dti-tk
- Fiber tracking: FACT from diffusion toolbox from trackvis
- Bundles segmentation: recoBundles from dipy
- Analysis: Statistical inference (VBA, tract-based statistics): dipy

## Installation instructions

### Linux/Mac
1 - Download and install Anaconda:
https://docs.anaconda.com/anaconda/install/linux/

2 - Navigate to ```./```
```
cd ./
```

3 - Create conda environment with the conda_environment.yml
file with all the necessary packages:
```
conda env create -f conda_environment.yml
```

IF there are any errors with cached-property, try installing first the conda environment only with Python. Then populate it with:

```
conda env update -f ./conda_environment.yml
```

4 - Activate the new conda environment:
```
conda activate conda_psycho_1_dwi_py367
```

4.1 - (Optional) In case of any update of the conda environment, run:
```
conda env update -f ./conda_environment.yml
```

5 - DSI Studio (PROBABLY it will only run in WINDOWS so don't follow these instructions for Ubuntu ---------

Install singularity:

https://github.com/hpcng/singularity/blob/master/INSTALL.md


To accept x11 GUI connections from docker images, run in the host sytem:
xhost +

Then run the docker image:
docker run -it --network=host --env DISPLAY=$DISPLAY  --env="DISPLAY=$DISPLAY" --env="QT_X11_NO_MITSHM=1" --privileged --volume="$HOME/.Xauthority:/root/.Xauthority:rw" -v /tmp/.X11-unix:/tmp/.X11-unix  -v /home/mario/workspace/tractography_dwi_analysis_cinac_david_mata/data/derivatives/qsiprep:/data --rm dsistudio/dsistudio:latest

(Finally, inside the docker image, run:)
dsi_studio

6 - fsl (optional):
https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/FslInstallation/ShellSetup

## Automatic configurations that can be manually tuned

---Configure heudiconv_subjs_anonymized_map.csv -----
This is a CSV in data/ which must have 2 columns (name, id) to assign an id to every participant.

---Configure heudiconv-----
Find what is the protocol name for every dicom file series. You can use DicomBrowser for this.
Then create "rules" in the heudiconv_config.py file to assign each DICOM series to a folder (wich will have the nifti file resultant of the conversion of dicom images).

If SDC with spin echo fieldmaps wants to be executed, the fieldmap PA and AP json sidecar must be edited to include:

"IntendedFor": [
    "ses-0/dwi/sub-{}_ses-0_run-01_dwi.nii.gz",
],

## Running instructions (preprocessing)

----- heudiconv -------

Outside the data/ folder run:
./heudiconv_command.sh

----- QSIPREP (q-space diffussion spectrum preprocessing) -------
Outside the data/ folder run:

qsiprep-docker data/bids data/derivatives -w data/work/ --fs-license-file /usr/local/freesurfer/license.txt \
--eddy-config eddy_params.json --output-resolution 1.2

docker run -ti --rm \
    -v /home/mario/workspace/tractography_dwi_analysis_cinac_david_mata/data/bids:/data:ro \
    -v /home/mario/workspace/tractography_dwi_analysis_cinac_david_mata/derivatives:/out \
    -v /usr/local/freesurfer:/usr/local/freesurfer \
    -v /home/mario/workspace/tractography_dwi_analysis_cinac_david_mata/:/eddy/ \
    pennbbl/qsiprep:latest \
    /data /out participant \
    --eddy-config /eddy/eddy_params.json --fs-license-file /usr/local/freesurfer/license.txt \
--output-resolution 1.2  -w data/work/

Consider using --dwi-only --ignore fieldmaps for faster tests (but after the tests, REMOVE IT)

CUDA Support is not implemented for RTX cards.

As of version 0.6.7 CUDA version 9.1 is supported in the QSIPrep container! To run locally using docker you will need the nvidia container runtime installed for Docker version 19.0.3 or higher. Singularity images will run with CUDA 9.1 with the -nv flag.

## Running instructions (processing)
Run the analysis/dti_pipeline.py file with the desired flags in the __main__ function at the end of the script.
