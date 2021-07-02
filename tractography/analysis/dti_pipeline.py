import subprocess
from pathlib import Path
from shutil import copyfile

import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
import scipy
import random
from matplotlib.ticker import FormatStrFormatter

import dipy.stats.analysis as dsa
import dipy.tracking.streamline as dts
from dipy.segment.clustering import QuickBundles
from dipy.segment.metric import (AveragePointwiseEuclideanMetric,
                                 ResampleFeature)
from dipy.data.fetcher import get_two_hcp842_bundles
import dipy.data as dpd
from dipy.io.streamline import load_trk
from dipy.io.image import load_nifti
from dipy.viz import window, actor
from dipy.data import fetch_bundles_2_subjects, read_bundles_2_subjects
from dipy.tracking.streamline import transform_streamlines
from dipy.io.image import load_nifti, save_nifti
from dipy.io.gradients import read_bvals_bvecs
from dipy.core.gradients import gradient_table
import dipy.reconst.dti as dti
from dipy.data import get_fnames
from dipy.segment.mask import median_otsu
from dipy.reconst.dti import (TensorModel, color_fa, fractional_anisotropy,
                              geodesic_anisotropy, mean_diffusivity,
                              axial_diffusivity, radial_diffusivity,
                              lower_triangular, mode as get_mode)
from dipy.io.streamline import load_trk, save_trk
from dipy.io.utils import create_tractogram_header
from dipy.io.stateful_tractogram import Space, StatefulTractogram
from dipy.viz import window, actor
from time import sleep
from dipy.data import two_cingulum_bundles
from dipy.align.streamlinear import StreamlineLinearRegistration
from dipy.tracking.streamline import set_number_of_points
from mne.stats import (ttest_1samp_no_p, bonferroni_correction, fdr_correction,
                       permutation_cluster_test, permutation_t_test, permutation_cluster_1samp_test)
from permutation import permutation_test
from statsmodels.stats.multitest import multipletests
import cluster
import corrstats
from sklearn.preprocessing import StandardScaler
import sklearn

import nibabel as nib
from nibabel import Nifti1Image
##########################################################
##########################################################
##########################################################

derivatives_data = os.path.join("..", "data", "derivatives")
dwi_data = os.path.join(derivatives_data, "qsiprep")
dti_data = os.path.join(derivatives_data, "dti")
tracts_data = os.path.join(derivatives_data, "tracts")
analysis_results_dir = "results"
all_results_dir = os.path.join(analysis_results_dir, "all")
all_results_tracts_dir = os.path.join(all_results_dir, "tracts")
all_results_vba_dir = os.path.join(all_results_dir, "vba")
all_results_corr_features_dir = os.path.join(all_results_tracts_dir, "correlations_features")
results_vba_dir = os.path.join(analysis_results_dir, "vba")
results_tracts_dir = os.path.join(analysis_results_dir, "tracts")
generic_afq_dir = os.path.join("tracts", "afq")
results_afq_dir = os.path.join(analysis_results_dir, generic_afq_dir)
generic_shape_dir = os.path.join("tracts", "shape")
results_shape_dir = os.path.join(analysis_results_dir, generic_shape_dir)

paths_proj = [dti_data, tracts_data, analysis_results_dir, results_vba_dir,
              results_tracts_dir, results_afq_dir, results_shape_dir,
              all_results_dir, all_results_vba_dir, all_results_tracts_dir,
              all_results_corr_features_dir]

for dir in paths_proj:
    if not os.path.exists(dir):
        os.mkdir(dir)

num_subjs = 0 # 0 = ALL
subjs_names = sorted([f.name for f in Path(dwi_data).iterdir() if f.is_dir()])
if num_subjs != 0:
    subjs_names = subjs_names[0:num_subjs]

def get_data(subj, name, type):
    path_dir = ""
    path = ""

    if type == "dwi":
        path_dir = os.path.join(dwi_data, subj, "ses-0", "dwi")
    elif type == "dti":
        path_dir = os.path.join(dti_data, subj)
    elif type == "tracts":
        path_dir = os.path.join(tracts_data, subj)

    if name == "dwi_dir" or name == "dti_dir":
        return path_dir

    if name == "dwi":
        path= "{}_ses-0_run-01_space-T1w_desc-preproc_dwi.nii.gz".format(subj)
    elif name == "bval":
        path= "{}_ses-0_run-01_space-T1w_desc-preproc_dwi.bval".format(subj)
    elif name == "bvec":
        path= "{}_ses-0_run-01_space-T1w_desc-preproc_dwi.bvec".format(subj)
    elif name == "mask":
        path= "{}_ses-0_run-01_space-T1w_desc-brain_mask.nii.gz".format(subj)
    elif name == "mask_mni":
        path = "{}_to_mni_tensor_mask.nii.gz".format(subj)
    elif name == "mask_mni_template":
        path = "{}_mask.nii.gz".format(subj)
    elif name == "tensor":
        path = "{}_tensor.nii.gz".format(subj)
    elif name == "tensor_mni":
        path = "{}_to_mni_tensor.nii.gz".format(subj)
    elif name == "fa":
        path = "{}_fa.nii.gz".format(subj)
    elif name == "fa_mni":
        path = "{}_to_mni_tensor_fa.nii.gz".format(subj)
    elif name == "tr_mni":
        path = "{}_to_mni_tensor_tr.nii.gz".format(subj)
    elif name == "whole_tract":
        path = "{}_whole_trk.trk".format(subj)
    elif name == "whole_tract_mni":
        path = "{}_whole_trk_mni.trk".format(subj)
    elif name == "rec_bundles_dir":
        path = "rec_bundles"
    elif name == "org_bundles_dir":
        path = "org_bundles"
    elif name == "rec_bundles":
        path = os.path.join("rec_bundles", "*.trk")
    elif name.startswith("tract__"):
        tract_name = name.split("__")[1]
        tract_filename = "{}_whole_trk_mni_{}__labels__recognized_orig.trk".format(subj, tract_name)
        path = os.path.join("org_bundles", tract_filename)
    elif name.startswith("tract_stats__"):
        tract_name = name.split("__")[1]
        tract_filename = "{}_whole_trk_mni_{}__labels__recognized_orig.trk.stat.txt".format(subj, tract_name)
        path = os.path.join("org_bundles", tract_filename)
    elif name.startswith("atlas_tract__"):
        tract_name = name.split("__")[1]
        tract_filename = "{}.trk".format(tract_name)
        path = os.path.join("rec_bundles", tract_filename)

    path = os.path.join(path_dir, path)

    return path

def get_output_dir(data_name):
    output_dir = ""

    if data_name == "afq":
        output_dir = generic_afq_dir
    elif data_name == "shape":
        output_dir = generic_shape_dir

    return output_dir

def correct_q_form(tensor_path, dti_original_path):
    tensor = nib.load(dti_original_path)
    tensor.set_qform(tensor.get_qform(), 1)
    nib.save(tensor, tensor_path)

    print("q_form correction applied")

    return 0

##########################################################
# ----------------- DWI preprocessing (denoise, correct Eddy currents, etc): qsiprep -----------------
##########################################################

# ----- heudiconv -------
#
# Outside the data/ folder run:
# ./heudiconv_command.sh
#
# ----- QSIPREP (q-space diffussion spectrum preprocessing) -------
# Outside the data/ folder run:
#
#
# docker run -ti --rm \
#     -v /home/mario/workspace/tractography_dwi_analysis_cinac_david_mata/data/bids:/data:ro \
#     -v /home/mario/workspace/tractography_dwi_analysis_cinac_david_mata/derivatives:/out \
#     -v /usr/local/freesurfer:/usr/local/freesurfer \
#     -v /home/mario/workspace/tractography_dwi_analysis_cinac_david_mata/:/eddy/ \
#     pennbbl/qsiprep:latest \
#     /data /out participant \
#     --eddy-config /eddy/eddy_params.json --fs-license-file /usr/local/freesurfer/license.txt \
# --output-resolution 1.2  -w data/work/ --ignore fieldmaps


##########################################################
# ----------------- DTI estimation (with seeds, end conditions, forbidden pahts, etc): dipy NLLS -----------------
##########################################################

# TODO: CSA: https://github.com/dipy/dipy/blob/c846bf5a23b7e95343a9cf231tract_df2653473602456/doc/interfaces/basic_flow.rst

def dti_estimation_nlls(fit_method="NLLS", overwrite=True):
    # https://github.com/dipy/dipy/blob/c846bf5a23b7e95343a9cf231tract_df2653473602456/dipy/workflows/reconst.py#L218

    global i, subj, dwi_path, bvec_path, mask_path, output_dir, p, line

    for i, subj in enumerate(subjs_names):
        print("----- Estimation DTI tensors for subj {} ({} of {})".format(subj, i + 1, len(subjs_names)))
        dwi_path = get_data(subj, "dwi", type="dwi")
        bvec_path = get_data(subj, "bvec", type="dwi")
        bval_path = get_data(subj, "bval", type="dwi")
        mask_path = get_data(subj, "mask", type="dwi")
        dti_output_dir = get_data(subj, "dti", type="dti")
        dti_output_path = get_data(subj, "tensor", type="dti")
        fa_output_path = get_data(subj, "fa", type="dti")

        if not overwrite and os.path.exists(dti_output_path):
            print("Skipping...already done")
            continue

        print("Loading input data...")
        data, affine = load_nifti(dwi_path)
        bvals, bvecs = read_bvals_bvecs(bval_path, bvec_path)
        gtab = gradient_table(bvals, bvecs)
        maskdata = load_nifti_data(mask_path).astype(bool)

        print("Fitting DTI with {}...".format(fit_method))
        tenmodel = dti.TensorModel(gtab, fit_method=fit_method)  # WLS NLLS
        tenfit = tenmodel.fit(data, maskdata)

        print("Saving DTI...")
        tensor_vals = lower_triangular(tenfit.quadratic_form)
        ten_img = nifti1_symmat(tensor_vals, affine=affine)
        if not os.path.exists(dti_output_dir):
            os.mkdir(dti_output_dir)
        nib.save(ten_img, dti_output_path)

        FA = fractional_anisotropy(tenfit.evals)
        FA[np.isnan(FA)] = 0
        FA = np.clip(FA, 0, 1)
        save_nifti(fa_output_path, FA.astype(np.float32), affine)

def load_nifti_data(fname, as_ndarray=True):
    """Load only the data array from a nifti file.
    Parameters
    ----------
    fname : str
        Full path to the file.
    as_ndarray: bool, optional
        convert nibabel ArrayProxy to a numpy.ndarray.
        If you want to save memory and delay this casting, just turn this
        option to False (default: True)
    Returns
    -------
    data: np.ndarray or nib.ArrayProxy
    See also
    --------
    load_nifti
    """
    img = nib.load(fname)
    return np.asanyarray(img.dataobj) if as_ndarray else img.dataobj

def make5d(input):
    """reshapes the input to have 5 dimensions, adds extra dimensions just
    before the last dimession
    """
    input = np.asarray(input)
    if input.ndim > 5:
        raise ValueError("input is already more than 5d")
    shape = input.shape
    shape = shape[:-1] + (1,)*(5-len(shape)) + shape[-1:]
    return input.reshape(shape)

def nifti1_symmat(image_data, *args, **kwargs):
    """Returns a Nifti1Image with a symmetric matrix intent
    Parameters
    -----------
    image_data : array-like
        should have lower triangular elements of a symmetric matrix along the
        last dimension
    all other arguments and keywords are passed to Nifti1Image
    Returns
    --------
    image : Nifti1Image
        5d, extra dimensions addes before the last. Has symmetric matrix intent
        code
    """
    image_data = make5d(image_data)
    last_dim = image_data.shape[-1]
    n = (np.sqrt(1+8*last_dim) - 1)/2
    if (n % 1) != 0:
        raise ValueError("input_data does not seem to have matrix elements")

    image = Nifti1Image(image_data, *args, **kwargs)
    hdr = image.header
    hdr.set_intent('symmetric matrix', (n,))
    return image

def dti_estimation(overwrite=False):
    global i, subj, dwi_path, bvec_path, mask_path, output_dir, p, line

    for i, subj in enumerate(subjs_names):
        print("----- Estimation DTI tensors for subj {} ({} of {})".format(subj, i + 1, len(subjs_names)))
        type_file = "dwi"
        dwi_path = get_data(subj, "dwi", type=type_file)
        bvec_path = get_data(subj, "bvec", type=type_file)
        bval_path = get_data(subj, "bval", type=type_file)
        mask_path = get_data(subj, "mask", type=type_file)
        output_dir = os.path.join(dti_data, subj)

        dipy_cmd = "dipy_fit_dti"
        dipy_output_dir = "--out_dir {}".format(output_dir)
        dipy_output_file= "--out_tensor {}_tensor.nii.gz".format(subj)
        dipy_extra_args = "--nifti_tensor"
        if overwrite:
            dipy_extra_args += " --force"

        dipy_full_cmd = "{} {} {} {} {} {} {} {}".format(
            dipy_cmd, dwi_path, bval_path, bvec_path, mask_path, dipy_output_dir, dipy_output_file, dipy_extra_args).\
            split(" ")

        run_cmd(dipy_full_cmd)

    print("DTI estimation done!")


def run_cmd(dipy_full_cmd):
    global p, line
    p = subprocess.Popen(dipy_full_cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)

    while True:
        line = p.stdout.readline()
        print(">>>", line)
        if not line: break


##########################################################
# ----------------- DTI postprocessing (Tensor normalization, TBSS, etc): tensor normalization with dti-tk -----------------
##########################################################

def tensor_normalization(preprocessing=True, force_population_mean=True, overwrite=False):
    global i, subj, subj_path

    tensor_names = [Path(get_data(subj, "tensor", type="dti")).name for subj in subjs_names]
    subjs_txt_filename = "subjs.txt"
    subjs_txt_path = os.path.join(dti_data, subjs_txt_filename)

    with open(subjs_txt_path, "w") as text_file:
        for tensor_name in tensor_names:
            text_file.write(tensor_name + "\n")

    if preprocessing:
        for i, subj in enumerate(subjs_names):
            print("----- Preprocessing for tensor normalization for subj {} ({} of {})".format(subj, i + 1, len(subjs_names)))
            subj_path = os.path.join(dwi_data, subj)
            subj_dti_dir = os.path.join(dti_data, subj)
            subj_dti_path = get_data(subj, "tensor", type="dti")

            dti_original_filename = "{}_tensor_original.nii.gz".format(subj)
            dti_original_path = os.path.join(subj_dti_dir, dti_original_filename)
            if not os.path.exists(dti_original_path):
                copyfile(subj_dti_path, dti_original_path)

            correct_q_form(subj_dti_path, dti_original_path)
            preprocessing_dti_tk(subj_dti_path, subj)

    population_mean_dir = os.path.join(dti_data, "population_mean")
    population_mean_dti_path = get_data("population_mean", "tensor", type="dti")
    #mean_subjs_path = os.path.join(population_mean_dir, "population_mean_tensor.nii.gz")
    if force_population_mean:
        population_mean_tensor(population_mean_dti_path, population_mean_dir, subjs_txt_path)

    for i, subj in enumerate(subjs_names):
        print("----- Tensor normalization for subj {} ({} of {})".format(subj, i + 1, len(subjs_names)))
        subj_path = os.path.join(dwi_data, subj)
        subj_dti_dir = os.path.join(dti_data, subj)
        population_mean_dir = os.path.join(dti_data, "population_mean")
        subj_dti_path = get_data(subj, "tensor", type="dti")
        population_mean_dti_path = get_data("population_mean", "tensor", type="dti")

        subj_to_population_tract_df_path = map_subj_space_to_mean_population_space(subj_dti_path, subj, subj_dti_dir,
                                                                             population_mean_dti_path, overwrite=overwrite)
        population_to_mni_tract_df_path = map_population_space_to_mni(population_mean_dti_path, "population_mean",
                                                                population_mean_dir, population_mean_dti_path, overwrite=overwrite)
        subj_to_mni_path = map_subj_space_to_mni(subj_dti_path, subj, subj_dti_dir, subj_to_population_tract_df_path,
                                                 population_to_mni_tract_df_path, overwrite=overwrite)
        extract_dti_scalars(subj_to_mni_path, subj_dti_dir, smooth=4, overwrite=overwrite)

    print("Tensor normalization done!")

def extract_dti_scalars(subj_to_mni_path, subj_dti_dir, smooth=4, overwrite=False):
    # --- Extract scalars ---
    subj_to_mni_name = Path(subj_to_mni_path).stem.split(".")[0]
    subj_to_mni_path_no_extension = os.path.join(subj_dti_dir, subj_to_mni_name)

    fa_path = "{}_fa.nii.gz".format(subj_to_mni_path_no_extension)
    if not overwrite and os.path.exists(fa_path):
        print("Skipping extract_dti_scalars...Already done")
        return fa_path

    dti_tk_cmd = "TVtool -in {} -fa -out {}".format(subj_to_mni_path, fa_path).split(" ")
    run_cmd(dti_tk_cmd)

    # TR (trace) = MD*3
    dti_tk_cmd = "TVtool -in {} -tr -out {}_tr.nii.gz".format(subj_to_mni_path, subj_to_mni_path_no_extension).\
        split(" ")
    run_cmd(dti_tk_cmd)

    dti_tk_cmd = "TVtool -in {} -ad -out {}_ad.nii.gz".format(subj_to_mni_path, subj_to_mni_path_no_extension).\
        split(" ")
    run_cmd(dti_tk_cmd)

    dti_tk_cmd = "TVtool -in {} -rd -out {}_rd.nii.gz".format(subj_to_mni_path, subj_to_mni_path_no_extension).\
        split(" ")
    run_cmd(dti_tk_cmd)

    # echo ~Smooth diffusion maps~
    dti_tk_cmd = "SVGaussianSmoothing -in {subj_path}_fa.nii.gz -fwhm {smooth} {smooth} {smooth} " \
                 "-out {subj_path}_fa{smooth}mm.nii.gz".format(
        subj_path=subj_to_mni_path_no_extension, smooth=smooth).split(" ")
    run_cmd(dti_tk_cmd)

    dti_tk_cmd = "SVGaussianSmoothing -in {subj_path}_tr.nii.gz -fwhm {smooth} {smooth} {smooth} " \
                 "-out {subj_path}_tr{smooth}mm.nii.gz".format(
        subj_path=subj_to_mni_path_no_extension, smooth=smooth).split(" ")
    run_cmd(dti_tk_cmd)

    dti_tk_cmd = "SVGaussianSmoothing -in {subj_path}_ad.nii.gz -fwhm {smooth} {smooth} {smooth} " \
                 "-out {subj_path}_ad{smooth}mm.nii.gz".format(
        subj_path=subj_to_mni_path_no_extension, smooth=smooth).split(" ")
    run_cmd(dti_tk_cmd)

    dti_tk_cmd = "SVGaussianSmoothing -in {subj_path}_rd.nii.gz -fwhm {smooth} {smooth} {smooth} " \
                 "-out {subj_path}_rd{smooth}mm.nii.gz".format(
        subj_path=subj_to_mni_path_no_extension, smooth=smooth).split(" ")
    run_cmd(dti_tk_cmd)

    print("FA maps and others done.")

def preprocessing_dti_tk(dti_path, subj):
    # ------- Correct diffusivity magnitude -------
    subj_dti_dir = os.path.join(dti_data, subj)

    dti_tk_cmd = "TVtool -in {dti_path} -scale 1000 -out {dti_path}".format(dti_path=dti_path).split(" ")
    run_cmd(dti_tk_cmd)
    # ------- Pre-processing before registration - ------
    dti_tk_cmd = "TVtool -in {} -norm".format(dti_path).split(" ")
    run_cmd(dti_tk_cmd)
    dti_norm_path = os.path.join(subj_dti_dir, "{}_tensor_norm.nii.gz".format(subj))
    dti_tk_cmd = "SVtool -in {} -stats".format(dti_norm_path).split(" ")
    run_cmd(dti_tk_cmd)
    # Check Whether DTI Volumes Share one Common Voxel Space
    dti_tk_cmd = "TVResample -in {} -align center -size 256 256 256 -vsize 1 1 1".format(dti_path).split(" ")
    run_cmd(dti_tk_cmd)
    dti_tk_cmd = "TVAdjustVoxelspace -in {} -origin 0 0 0".format(dti_path).split(" ")
    run_cmd(dti_tk_cmd)


def map_subj_space_to_mni(dti_path, subj, dti_dir, subj_to_population_tract_df_path, population_to_mni_tract_df_path, overwrite=False):
    subj_mni_path = get_data(subj, "tensor_mni", type="dti")
    dti_name = Path(subj_mni_path).stem.split(".")[0]
    subj_mni_tract_df_path = dti_name + ".tract_df.nii.gz"
    subj_mni_tract_df_path = os.path.join(dti_dir, subj_mni_tract_df_path)
    subj_mni_mask_path = get_data(subj, "mask_mni", type="dti")
    mni_template_path = os.path.join(dti_data, "mni_template", "mni_template.nii.gz")
    subj_mni_tr_name = Path(subj_mni_path).stem.split(".")[0] + "_tr.nii.gz"
    subj_mni_tr_path = os.path.join(dti_dir, subj_mni_tr_name)

    if not overwrite and os.path.exists(subj_mni_path):
        print("Skipping map_subj_space_to_mni...Already done")
        return subj_mni_path

    # ------- Mapping from the subj space to the standard (MNI) space -------
    dti_tk_cmd = "tract_dfComposition -tract_df2 {} -tract_df1 {} -out {}".format(
        subj_to_population_tract_df_path, population_to_mni_tract_df_path, subj_mni_tract_df_path).split(" ")
    run_cmd(dti_tk_cmd)

    # ------- Applying the mapping from the subj space to the standard (MNI) space -------
    dti_tk_cmd = "deformationSymTensor3DVolume -in {} -trans {} -target {} -out {}".format(
        dti_path, subj_mni_tract_df_path, mni_template_path, subj_mni_path).split(" ")
    run_cmd(dti_tk_cmd)

    dti_tk_cmd = "TVtool -tr -in {}".format(subj_mni_path).split(" ")
    run_cmd(dti_tk_cmd)

    dti_tk_cmd = "BinaryThresholdImageFilter {} {} 0 .01 0 1 0".format(
        subj_mni_tr_path, subj_mni_mask_path).split(" ")
    run_cmd(dti_tk_cmd)

    print("Mapping done! The output file {} is now the subj dti data registered into the "
          "MNI template data".format(subj_mni_path))

    return subj_mni_path

def map_population_space_to_mni(dti_path, subj, dti_dir, population_mean_dti_path, overwrite=False):
    # Download MNI atlas and its mask:
    # https://www.nitrc.org/frs/download.php/11301/IITmean_tensor_256.nii.gz
    # https://www.nitrc.org/frs/download.php/11302/IITmean_tensor_mask_256.nii.gz
    dti_name = Path(dti_path).stem.split(".")[0]
    transf_name = dti_name + "_aff.nii.gz"
    transf_path = os.path.join(dti_dir, transf_name)
    mni_template_path =  os.path.join(dti_data, "mni_template", "mni_template.nii.gz")
    mask_name = dti_name + "_mask.nii.gz"
    mean_mask_path = os.path.join(dti_dir, mask_name)

    transf_filename = Path(transf_path).stem.split(".")[0]
    diffeo_name = transf_filename + "_diffeo.tract_df.nii.gz"
    diffeo_path = os.path.join(dti_dir, diffeo_name)
    population_to_mni_path = os.path.join(dti_dir, "{}_to_mni_tensor.tract_df.nii.gz".format(subj))
    transf_aff_name = dti_name + ".aff"
    transf_aff_path = os.path.join(dti_dir, transf_aff_name)

    if not overwrite and os.path.exists(population_to_mni_path):
        print("Skipping map_population_space_to_mni...Already done")
        return population_to_mni_path

    dti_tk_cmd = "dti_rigid_reg {} {} EDS 4 4 4 0.001".format(
        mni_template_path, population_mean_dti_path).split(" ")
    run_cmd(dti_tk_cmd)

    dti_tk_cmd = "dti_affine_reg {} {} EDS 4 4 4 0.001 1".format(
        transf_path, population_mean_dti_path).split(" ")
    run_cmd(dti_tk_cmd)

    dti_tk_cmd = "TVtool -tr -in {}".format(transf_path).split(" ")
    run_cmd(dti_tk_cmd)

    transf_tr_name = Path(transf_path).stem.split(".")[0] + "_tr.nii.gz"
    transf_tr_path = os.path.join(dti_dir, transf_tr_name)
    dti_tk_cmd = "BinaryThresholdImageFilter {} {} 0 .01 100 1 0".format(
        transf_tr_path, mean_mask_path).split(" ")
    run_cmd(dti_tk_cmd)

    iters = 1
    dti_tk_cmd = "dti_diffeomorphic_reg {} {} {} 1 {} 0.002".format(
        mni_template_path, transf_path, mean_mask_path, iters).split(" ")
    run_cmd(dti_tk_cmd)

    dti_tk_cmd = "tract_dfRightComposeAffine -aff {} -tract_df {} -out {}".format(
        transf_aff_path, diffeo_path, population_to_mni_path).split(" ")
    run_cmd(dti_tk_cmd)

    print("Registration done! The output file {} is now the registration of"
          " the population mean dti data into the MNI data. IT IS NOT the mapped data, it is only the transformation."
          .format(population_to_mni_path))

    return population_to_mni_path

def map_subj_space_to_mean_population_space(subj_dti_path, subj, subj_dti_dir, population_mean_dti_path, overwrite=False):
    # Step 1: Create a subjs.txt file.
    subj_dti_name = Path(subj_dti_path).stem.split(".")[0]

    mask_name = subj_dti_name + "_mask.nii.gz"
    mask_path = os.path.join(subj_dti_dir, mask_name)

    population_mean_dir = os.path.join(dti_data, "population_mean")
    subj_to_population_tract_df_path = os.path.join(subj_dti_dir, "{}_to_population.tract_df.nii.gz".format(subj))
    transf_name = subj_dti_name + "_aff.nii.gz"
    transf_path = os.path.join(subj_dti_dir, transf_name)
    transf_aff_name = subj_dti_name + ".aff"
    transf_aff_path = os.path.join(subj_dti_dir, transf_aff_name)

    if not overwrite and os.path.exists(subj_to_population_tract_df_path):
        print("Skipping map_subj_space_to_mean_population_space...Already done")
        return subj_to_population_tract_df_path

    # Step 2: Mean instead of an already existing template_popoulation (DONE outside of this function)

    # Step 3: Rigid Alignment of DTI Volumes
    dti_tk_cmd = "dti_rigid_reg {} {} EDS 4 4 4 0.01".format(population_mean_dti_path, subj_dti_path).split(" ")
    run_cmd(dti_tk_cmd)

    # Step 4: Affine Alignment of DTI Volumes
    dti_tk_cmd = "dti_affine_reg {} {} EDS 4 4 4 0.01 1".format(transf_path, subj_dti_path).split(" ")
    run_cmd(dti_tk_cmd)

    # Step 5a: Create subj dti mask:
    dti_tk_cmd = "TVtool -tr -in {}".format(transf_path).split(" ")
    run_cmd(dti_tk_cmd)
    transf_tr_name = Path(transf_path).stem.split(".")[0] + "_tr.nii.gz"
    transf_tr_path = os.path.join(subj_dti_dir, transf_tr_name)
    dti_tk_cmd = "BinaryThresholdImageFilter {} {} 0 .01 100 1 0".format(transf_tr_path, mask_path).split(" ")
    run_cmd(dti_tk_cmd)

    # Step 5b: Deformable Alignment of DTI Volumes
    iters = 1
    dti_tk_cmd = "dti_diffeomorphic_reg {} {} {} 1 {} 0.002".format(
        population_mean_dti_path, transf_path, mask_path, iters).split(" ")
    run_cmd(dti_tk_cmd)

    # Step 6: Combined displacement field from the native space to the template space
    diffeo_name = Path(transf_path).stem.split(".")[0] + "_diffeo.tract_df.nii.gz"
    diffeo_path = os.path.join(subj_dti_dir, diffeo_name)
    dti_tk_cmd = "tract_dfRightComposeAffine -aff {} -tract_df {} -out {}".format(
        transf_aff_path, diffeo_path, subj_to_population_tract_df_path).split(" ")
    run_cmd(dti_tk_cmd)

    # Step 7: Mapping the subject data to the template space (Not necessary since then we map subj to MNI using this population transformation)
    # dti_tk_cmd = "deformationSymTensor3DVolume -in ${subj}_tensor.nii.gz -trans ${subj}_to_population.tract_df.nii.gz -target population_mean_tensor.nii.gz -out ${subj}_to_mean.nii.gz"
    # run_cmd(dti_tk_cmd)

    print("Mapping done! The output file {} is now the subj dti data registered into"
          " the population template data".format(subj_to_population_tract_df_path))

    return subj_to_population_tract_df_path


def population_mean_tensor(population_mean_dti_path, population_mean_dir, subjs_txt_path):
    if not os.path.exists(population_mean_dir):
        os.mkdir(population_mean_dir)
    for i, this_subj in enumerate(subjs_names):
        this_dti_path = get_data(this_subj, "tensor", type="dti")
        this_dti_name = Path(this_dti_path).name
        copyfile(this_dti_path, this_dti_name)
    dti_tk_cmd = "TVMean -in {} -out {}".format(subjs_txt_path, population_mean_dti_path).split(" ")
    run_cmd(dti_tk_cmd)
    dti_tk_cmd = "TVResample -in {} -vsize 1 1 1 -size 256 256 256".format(population_mean_dti_path).split(" ")
    run_cmd(dti_tk_cmd)
    # Cleanup
    for i, this_subj in enumerate(subjs_names):
        this_dti_path = get_data(this_subj, "tensor", type="dti")
        this_dti_name = Path(this_dti_path).name
        os.remove(this_dti_name)


##########################################################
# ----------------- Fiber tracking: FACT from diffusion toolbox from trackvis  -----------------
##########################################################

def fiber_tracking():
    global overwrite, i, subj, subj_path, mask_path, output_dir, p, line
    overwrite = False

    diff_toolkit_path = "/home/mario/Downloads/Diffusion_Toolkit_v0.6.4.1_x86_64/dtk"
    for i, subj in enumerate(subjs_names):
        print("----- Fiber tracking for subj {} ({} of {})".format(subj, i + 1, len(subjs_names)))
        dti_dir = get_data(subj, "dti_dir", type="dti")
        dti_name = Path(get_data(subj, "tensor_mni", type="dti"))
        dti_name_no_extension = dti_name.stem.split("_tensor")[0]
        dti_path = os.path.join(dti_dir, dti_name_no_extension)
        mask_path = get_data(subj, "mask_mni", type="dti")
        #mask_path = get_data("mni_template", "mask_mni_template", type="dti")
        fa_path = get_data(subj, "fa_mni", type="dti")

        diff_toolk_cmd = os.path.join(diff_toolkit_path, "dti_tracker")
        diff_toolk_output = get_data(subj, "whole_tract", type="tracts")
        subj_tracts_path = os.path.join(tracts_data, subj)
        if not os.path.exists(subj_tracts_path):
            os.mkdir(subj_tracts_path)
        diff_toolk_mask = "-m {}".format(mask_path)  # ${subj}_tensor_mask.nii.gz"
        diff_toolk_mask_2 = "-m2 {}".format(fa_path)  # ${subj}_tensor_mask.nii.gz"
        diff_toolk_filetype = "-it nii.gz"
        diff_toolk_algorithm = "-fact" # "-fact"

        diff_toolk_full_cmd = "{} {} {} {} {} {}".format(
            diff_toolk_cmd, dti_path, diff_toolk_output, diff_toolk_mask,
            diff_toolk_filetype, diff_toolk_algorithm).split(" ")

        run_cmd(diff_toolk_full_cmd)

    print("Fiber tracking done!")

##########################################################
# ----------------- Bundles segmentation: recoBundles from dipy -----------------
##########################################################

# https://dipy.org/documentation/1.4.0./interfaces/bundle_segmentation_flow/

# Download from: https://ndownloader.figshare.com/files/26842853

# dipy_horizon "subj0_to_mni_trck.trk" "atlas.trk" --random_color

# dipy_slr "atlas.trk" "subj0_to_mni_trck.trk" --force

# dipy_recobundles "moved.trk" "atlas/bundles/*.trk" --force --mix_names --out_dir "rb_output"

# dipy_horizon moved_CCMid__recognized.trk ../atlas/bundles/CCMid.trk --random_colors

# Output of RecoBundles will be in native space. To get bundles in subject’s original space, run following commands:

# mkdir org_output

# dipy_labelsbundles 'subj0_to_mni_trck.trk' 'rb_output/*.npy' --mix_names --out_dir "org_output"

def slr(subj, num_points_per_streamline=20):
    # ---- SLR (Streamline-Based Linear Registration) -----
    print("-- Starting SLR (Streamline-Based Linear Registration)")
    global output_dir, p, line
    overwrite = True

    atlas_whole_trk = get_data("atlas", "whole_tract", type="tracts")
    subj_trck = get_data(subj, "whole_tract", type="tracts")
    subj_mni_trck = Path(get_data(subj, "whole_tract_mni", type="tracts")).name
    output_dir = os.path.join(tracts_data, subj)

    # TODO: All streamlines need to have the same number of points.
    # Use dipy.tracking.streamline.set_number_of_points to adjust your streamlines

    # #---------------------------------------------
    # subj_mni_tensor_path = get_data(subj, "tensor_mni", type="dti")
    # atlas_mask_mni_tensor_path = get_data("mni_template", "mask_mni", type="dti")
    # # subj_tract_path = get_data(subj, "whole_tract", type="tracts")
    # atlas_tract_path = get_data("atlas", "whole_tract", type="tracts")
    # # subj_tract = load_trk(subj_tract_path, "same", bbox_valid_check=False).streamlines
    # atlas_tract = load_trk(atlas_tract_path, "same", bbox_valid_check=False)
    #
    # # num_points_subj = [len(streamline) for streamline in subj_tract]
    # num_points_atlas = [len(streamline) for streamline in atlas_tract]
    #
    # # cb_subj = set_number_of_points(subj_tract, 20)
    # cb_atlas = set_number_of_points(atlas_tract.streamlines, 20)
    #
    # # sft_subj = StatefulTractogram(cb_subj, subj_mni_tensor_path, Space.RASMM)
    # # save_trk(sft_subj, "whole_tract_reoriented.trk", subj_tract)
    #
    # img = nib.Nifti1Image(np.zeros((2, 2, 2)), atlas_tract.affine)
    # # Move the streamlines to this other space, and report it:
    # sft_atlas = StatefulTractogram(cb_atlas, img, Space.RASMM)
    # save_trk(sft_atlas, "atlas_tract_reoriented.trk", atlas_tract)

    #---------------------------------------------

    dipy_cmd = "dipy_slr"
    dipy_output_dir = "--out_dir {}".format(output_dir)
    dipy_output_file = "--out_moved {}".format(subj_mni_trck)
    dipy_num_points_per_streamline = "--nb_pts {}".format(num_points_per_streamline)
    dipy_extra_args = ""
    if overwrite:
        dipy_extra_args += "--force"

    dipy_full_cmd = "{} {} {} {} {} {} {}".format(
        dipy_cmd, atlas_whole_trk, subj_trck, dipy_output_dir, dipy_output_file,
        dipy_num_points_per_streamline, dipy_extra_args). \
        split(" ")
    run_cmd(dipy_full_cmd)

    print("SLR done!")

def reco_bundles(subj):
    # https://dipy.org/documentation/1.4.0./examples_built/bundle_registration/#example-bundle-registration
    # https://dipy.org/documentation/1.1.0./reference_cmd/dipy_recobundles

    # ---- Bundles segmentation (dipy reboBundles) -----
    print("-- Starting Bundles segmentation (dipy reboBundles)")
    overwrite = True
    global output_dir, p, line

    subj_trck = get_data(subj, "whole_tract", type="tracts")
    subj_mni_trck = get_data(subj, "whole_tract_mni", type="tracts")
    atlas_bundles = get_data("atlas", "rec_bundles", type="tracts")
    subj_bundles_dir = get_data(subj, "rec_bundles_dir", type="tracts")
    subj_org_bundles_dir = get_data(subj, "org_bundles_dir", type="tracts")
    output_dir = os.path.join(dti_data, subj)

    dipy_cmd = "dipy_recobundles"
    dipy_output_dir = "--out_dir {}".format(subj_bundles_dir)
    dipy_extra_args = "--mix_names "
    if overwrite:
        dipy_extra_args += "--force"

    dipy_full_cmd = "{} {} {} {} {}".format(
        dipy_cmd, subj_mni_trck, atlas_bundles, dipy_output_dir, dipy_extra_args). \
        split(" ")
    run_cmd(dipy_full_cmd)

    # --------------------
    dipy_cmd = "dipy_labelsbundles"
    dipy_input_bundles_transf = "{}/*.npy".format(subj_bundles_dir)
    dipy_output_dir = "--out_dir {}".format(subj_org_bundles_dir)
    dipy_extra_args = "--mix_names "
    if overwrite:
        dipy_extra_args += "--force"
    dipy_full_cmd = "{} {} {} {} {}".format(
        dipy_cmd, subj_trck, dipy_input_bundles_transf, dipy_output_dir, dipy_extra_args). \
        split(" ")
    run_cmd(dipy_full_cmd)

    filenames = os.listdir(subj_org_bundles_dir)
    for old_filename in filenames:
        if old_filename.endswith("_orig.trk"):
            old_filepath = os.path.join(subj_org_bundles_dir, old_filename)
            try:
                new_filename = subj + "_whole_trk_" + old_filename.split("whole_trk_")[2]
                new_filepath = os.path.join(subj_org_bundles_dir, new_filename)
                os.rename(old_filepath, new_filepath)
            except IndexError as e:
                pass

    print("reco_bundles done!")

def bundles_segmentation(run_slr=True, run_reco_bundles=True):
    # https://dipy.org/documentation/1.4.0./reference/dipy.segment/#recobundles
    global i, subj, subj_path, output_dir, p, line

    for i, subj in enumerate(subjs_names):
        print("----- Bundles segmentation for subj {} ({} of {})".format(subj, i + 1, len(subjs_names)))
        subj_path = os.path.join(dwi_data, subj)

        # ---- SLR (Streamline-Based Linear Registration) -----
        if run_slr:
            slr(subj)

        # ---- Bundles segmentation (dipy reboBundles) -----
        if run_reco_bundles:
            reco_bundles(subj)

    print("Bundles segmentation done!")

##########################################################
# ----------------- Analysis: Statistical inference (VBA, tract-based statistics): dipy -----------------
##########################################################

def vba(merge_fa=True, smooth=4, num_permutations=1000):
    # http://www.diffusion-imaging.com/2015/10/dti-tutorial-2-normalization-and.html
    # http://brainimaging.waisman.wisc.edu/~tromp/example_data/DTI_Tutorial_2_answers.ptract_df
    # https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/GLM/CreatingDesignMatricesByHand

    design_matrix_txt_path = os.path.join(results_vba_dir, "design_matrix.txt")
    design_matrix_mat_path = os.path.join(results_vba_dir, "design_matrix.mat")
    contrasts_matrix_txt_path = os.path.join(results_vba_dir, "contrasts_matrix.txt")
    contrasts_matrix_con_path = os.path.join(results_vba_dir, "contrasts_matrix.con")
    tensor_mni_mask = get_data(subjs_names[0], "mask_mni", "dti")
    output_map_path = os.path.join(results_vba_dir, "randomise_out.nii.gz")

    print("------- Starting vba analysis with FSL randomise ---------")

    # echo ~~~Merge FA images~~~
    print("Running FSL merge...")
    all_merge_fa_name = "all_merge_tensor_mni_FA.nii.gz"
    all_merge_fa_path = os.path.join(dti_data, all_merge_fa_name)
    fa_paths = []
    for j, subj in enumerate(subjs_names):
        fa_path = get_data(subj, "fa_mni", type="dti")
        fa_paths.append(fa_path)

    fsl_cmd = "fslmerge"

    merge_mode = "-t"
    fa_paths_str = " ".join(fa_paths)
    input_fa = "{}".format(fa_paths_str)
    fsl_full_cmd = "{} {} {} {}".format(
        fsl_cmd, merge_mode, all_merge_fa_path, input_fa). \
        split(" ")
    if merge_fa:
        run_cmd(fsl_full_cmd)

    print("Converting design_matrix.txt and contrasts_matrix.txt to FSL format (.mat and .con)...")
    fsl_cmd = "Text2Vest"
    fsl_full_cmd = "{} {} {}".format(
        fsl_cmd, design_matrix_txt_path, design_matrix_mat_path). \
        split(" ")
    run_cmd(fsl_full_cmd)

    fsl_full_cmd = "{} {} {}".format(
        fsl_cmd, contrasts_matrix_txt_path, contrasts_matrix_con_path). \
        split(" ")
    run_cmd(fsl_full_cmd)

    print("Running FSL randomise...")
    fsl_cmd = "randomise"
    input_fa_merge = "-i {}".format(all_merge_fa_path)
    output_map = "-o {}".format(output_map_path)
    design_matrix = "-d {}".format(design_matrix_mat_path)
    contrasts_matrix = "-t {}".format(contrasts_matrix_con_path)
    smooth_fwhm = "-v {}".format(smooth)
    num_permutations = "-n {}".format(num_permutations)
    mask = "-m {}".format(tensor_mni_mask)
    algorithm = "-C 1.6"  # -T # -C

    fsl_full_cmd = "{} {} {} {} {} {} {} {} {}".format(
        fsl_cmd, input_fa_merge, output_map, design_matrix, contrasts_matrix, mask, smooth_fwhm, num_permutations, algorithm). \
        split(" ")
    run_cmd(fsl_full_cmd)


# ----------------- AFQ tract profiles -----------------
# https://dipy.org/documentation/1.4.0./examples_built/afq_tract_profiles/#example-afq-tract-profiles

def afq_profiles(show_interactive_bundle=False, measure="fa", overwrite=False):
    atlas_bundles_dir = get_data("atlas", "rec_bundles_dir", type="tracts")
    tracts_names = sorted([f.stem for f in Path(atlas_bundles_dir).iterdir() if f.suffix == ".trk"])

    for i, tract_name in enumerate(tracts_names):
        print("----- AFQ profile ({}) for {} ({} of {})".format(measure,
              tract_name, i + 1, len(tracts_names)))
        all_tract_stats = {}
        tract_df_name = "{}__{}.csv".format(tract_name, measure)
        tract_df_path = os.path.join(results_afq_dir, tract_df_name)
        for j, subj in enumerate(subjs_names):
            print("----- Generating profile tract for subj {} ({} of {})".format(subj, j + 1, len(subjs_names)))
            subj_bundles_dir = get_data(subj, "org_bundles_dir", type="tracts")
            subj_tract_path = get_data(subj, "tract__{}".format(tract_name), type="tracts")
            atlas_tract_path = get_data("atlas", "atlas_tract__{}".format(tract_name), type="tracts")

            try:
                afq_name = "{}_{}_profile.csv".format(tract_name, measure)
                afq_path = os.path.join(subj_bundles_dir, afq_name)
                if not overwrite and os.path.exists(afq_path):
                    print("Skipping...already done")
                    afq_profile_tract_df = pd.read_csv(afq_path, index_col=0)
                else:
                    subj_tract = load_trk(subj_tract_path, "same", bbox_valid_check=False).streamlines
                    atlas_tract = load_trk(atlas_tract_path, "same", bbox_valid_check=False).streamlines

                    # ===========================================================
                    feature = ResampleFeature(nb_points=100)
                    metric = AveragePointwiseEuclideanMetric(feature)
                    # ===========================================================
                    qb = QuickBundles(np.inf, metric=metric)
                    cluster_tract = qb.cluster(atlas_tract)
                    standard_tract = cluster_tract.centroids[0]
                    # ===========================================================
                    oriented_tract = dts.orient_by_streamline(subj_tract, standard_tract)
                    # ===========================================================
                    scalar_mni_path = get_data(subj, "{}_mni".format(measure), type="dti")
                    scalar, scalar_affine = load_nifti(scalar_mni_path)
                    # ===========================================================

                    w_tract = dsa.gaussian_weights(oriented_tract)  # oriented_tract
                    afq_profile = dsa.afq_profile(scalar, oriented_tract, scalar_affine,  # oriented_tract
                                                   weights=w_tract)
                    afq_profile_tract_df = pd.DataFrame(afq_profile, columns=[measure.upper()])
                    afq_profile_tract_df.to_csv(afq_path)

                    # Visualization

                    ylabel = measure
                    xlabel = "Segments along {} tract".format(tract_name)
                    fig_name = "{}_{}_profile.png".format(tract_name, measure)
                    fig_path = os.path.join(subj_bundles_dir, fig_name)
                    fig, ax = plt.subplots()
                    ax.set_title(subj)
                    ax.set_ylabel(ylabel)
                    ax.plot(afq_profile)
                    ax.set_xlabel(xlabel)
                    fig.savefig(fig_path)
                    plt.show()

                    show_bundle(scalar, scalar_affine, subj_tract, tract_name, subj_bundles_dir,
                                show_interactive_bundle=show_interactive_bundle)

                # tract_stats = {
                #     "1q_mean": afq_profile_tract_df.iloc[0:24, 0].mean(),
                #     "2q_mean": afq_profile_tract_df.iloc[25:50, 0].mean(),
                #     "3q_mean": afq_profile_tract_df.iloc[51:75, 0].mean(),
                #     "4q_mean": afq_profile_tract_df.iloc[76:99, 0].mean(),
                #     "global_mean": afq_profile_tract_df.mean(),
                # }
                # tract_stats_tract_df = pd.DataFrame(tract_stats)
                # all_tract_stats[subj] = tract_stats_tract_df
                all_tract_stats[subj] = afq_profile_tract_df.T
            except Exception as e:
                filename = "{}_error_afq.log".format(tract_name)
                error_msg = "Failed to generate AFQ profile for tract {} of subj {}. \n " \
                            "This error was probably caused by an incorrect bundle recognition \n {}".format(
                    tract_name, subj, e)
                print(error_msg)
                filepath = os.path.join(subj_bundles_dir, filename)
                with open(filepath, "w") as text_file:
                    text_file.write(error_msg)

        if not all_tract_stats:
            filename = "{}__{}_error.log".format(measure, tract_name)
            error_msg = "Failed to generate stats for AFQ profile for tract {}.  \n " \
                        "This tract was not recognized in any subject.".format(
                tract_name)
            print(error_msg)
            filepath = os.path.join(results_afq_dir, filename)
            with open(filepath, "w") as text_file:
                text_file.write(error_msg)
        else:
            all_tract_stats_tract_df = pd.concat(all_tract_stats, axis=0, ignore_index=True, sort=True).astype("float")
            all_tract_stats_tract_df.insert(0, "subj_id", all_tract_stats.keys())
            all_tract_stats_tract_df.insert(0, "group", [0]*all_tract_stats_tract_df.shape[0])
            all_tract_stats_tract_df.to_csv(tract_df_path, index=False)

def bundles_shape_analysis(dsi_studio_cmd_ran=False, overwrite=False):
    global i, subj, subj_path, output_dir, p, line

    if not dsi_studio_cmd_ran:
        print("...dsi_studio analysis command has not been executed yet. Follow these instructions:")
        print("======= Instructions to generate bundle shape stats ========")
        project_path = "$PWD" + os.sep + ".."
        dsi_studio_cmd = "docker run -ti --rm -e DISPLAY=$DISPLAY -v /tmp/.X11-unix:/tmp/.X11-unix -v $PWD/..:/data dsistudio/dsistudio:latest".format(project_path)
        print("From this analysis folder ({}), run:".format(os.getcwd()))
        print(dsi_studio_cmd)

        dsi_studio_cmd_py = "python3 data/analysis/dsi_studio_tracts_shape_analysis.py"
        print("Now inside the docker bash console, run:")
        print(dsi_studio_cmd_py)
        print("When it finishes, set to true the  dsi_studio_cmd_ran argument of this function")
        print("==========================================")

        return -1
    else:
        atlas_bundles_dir = get_data("atlas", "rec_bundles_dir", type="tracts")
        tracts_names = sorted([f.stem for f in Path(atlas_bundles_dir).iterdir() if f.suffix == ".trk"])
        all_tract_stats_tract_df = {}

        for i, tract_name in enumerate(tracts_names):
            print("----- bundles_shape_analysis {} ({} of {})".format(tract_name, i + 1, len(tracts_names)))
            all_tract_stats = {}
            tract_df_name = "{}__shape.csv".format(tract_name)
            tract_df_path = os.path.join(results_shape_dir, tract_df_name)

            if not overwrite and os.path.exists(tract_df_path):
                print("Skipping...already done")
                continue

            for j, subj in enumerate(subjs_names):
                print("----- bundles_shape_analysis for subj {} ({} of {})".format(subj, j + 1, len(subjs_names)))
                tract_stats_path = get_data(subj, "tract_stats__{}".format(tract_name), type="tracts")

                if os.path.exists(tract_stats_path):
                    tract_stats_tract_df = pd.read_csv(tract_stats_path, sep='\t', index_col=False, header=None).T
                    tract_stats_tract_df.columns = tract_stats_tract_df.iloc[0]
                    tract_stats_tract_df = tract_stats_tract_df.drop(tract_stats_tract_df.index[0])

                    all_tract_stats[subj] = tract_stats_tract_df

            if not all_tract_stats:
                filename = "{}__shape_error.log".format(tract_name)
                error_msg = "Failed to generate stats for shape profile for tract {}.  \n " \
                            "This tract was not recognized in any subject.".format(
                    tract_name)
                print(error_msg)
                filepath = os.path.join(results_shape_dir, filename)
                with open(filepath, "w") as text_file:
                    text_file.write(error_msg)
            else:
                all_tract_stats_tract_df = pd.concat(all_tract_stats, axis=0,  ignore_index=True, sort=True).astype("float")
                all_tract_stats_tract_df.insert(0, "subj_id", all_tract_stats.keys())
                all_tract_stats_tract_df.insert(0, "group", [0]*all_tract_stats_tract_df.shape[0])
                all_tract_stats_tract_df.to_csv(tract_df_path, index=False)

        return all_tract_stats_tract_df

def show_bundle(fa, fa_affine, subj_tract, tract_name, subj_bundles_dir, show_interactive_bundle=False):
    bundle_native = transform_streamlines(subj_tract, np.linalg.inv(fa_affine))
    scene = window.Scene()
    stream_actor = actor.line(bundle_native, fa, linewidth=0.1)
    bar = actor.scalar_bar()
    scene.add(stream_actor)
    scene.add(bar)

    if show_interactive_bundle:
        window.show(scene, size=(600, 600), reset_camera=False)

    fig_name = "{}.png".format(tract_name)
    fig_path = os.path.join(subj_bundles_dir, fig_name)
    window.record(scene, out_path=fig_path, size=(600, 600))

def stats_t_test(data_name="afq", only_hc_pd=False, test_type="permutation",
                 measure="fa", correct_after_perms=True, fontsize=16, overwrite=False):
    data_dir = get_output_dir(data_name)
    output_dir = os.path.join(all_results_dir, data_dir)
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)

    group_dirs = [f.name for f in Path(analysis_results_dir).iterdir() if f.is_dir() and f.name.startswith("group_")]

    all_subjs_results = {}
    for group_dir in group_dirs:
        subj_results_dir = os.path.join(analysis_results_dir, group_dir, data_dir)

        group = group_dir.split("_")[1]
        print("---- Reading {} group data".format(group))

        filenames = sorted([f.name for f in Path(subj_results_dir).iterdir() if f.suffix == ".csv" and
                            not "__stats_p_vals" in f.name and "__{}".format(measure) in f.name])
        for file_name in filenames:
            print("Reading data for {}".format(file_name))
            tract_df_path = os.path.join(subj_results_dir, file_name)
            tract_df = pd.read_csv(tract_df_path, index_col=False)
            file_name_no_ext = Path(file_name).stem

            tract_df["group"] = [group] * tract_df.shape[0]

            if file_name_no_ext not in all_subjs_results:
                all_subjs_results[file_name_no_ext] = {}

            if only_hc_pd and group != "hc":
                group = "pd"

            if group in all_subjs_results[file_name_no_ext]:
                tract_df_group = all_subjs_results[file_name_no_ext][group]
                all_subjs_results[file_name_no_ext][group] = pd.concat([tract_df_group, tract_df])
            else:
                all_subjs_results[file_name_no_ext][group] = tract_df

    for tract_name, subjs_vals in all_subjs_results.items():
        tract_df = pd.concat(subjs_vals, axis=0, ignore_index=True, sort=False)

        column_list = [x for x in tract_df.columns if x != 'group' and x != "subj_id"]
        group_names = tract_df["group"].unique().tolist()
        if len(group_names) == 2:
            group_names_alt = ["Hypersexual", "Non-hypersexual"]
            group_colors = ["#ff2600", "#73fdff"]  # Red and turquoise
        else:
            group_names_alt = ["HC"]
            group_colors = ["blue"]

        # --------------
        if data_name == "afq":
            if "fa" in measure:
                ylabel = "FA"
            elif "tr" in measure:
                ylabel = "MD*"
            else:
                ylabel = measure

            xlabel = "Segments along {} tract".format(tract_name)
            fig_name = "{}_profile_comparison.png".format(tract_name)
            fig_path = os.path.join(output_dir, fig_name)
            fig, ax = plt.subplots()
            # ax.set_title("All subjs")
            ax.set_ylabel(ylabel, fontsize=fontsize)
            # ax.set_xlabel(xlabel, fontsize=10)
            data_plot_mean = tract_df.groupby(['group']).mean()
            data_plot_std = tract_df.groupby(['group']).std()

        tract_already_done = False
        for num_group, group_i in enumerate(group_names):
            if data_name == "afq":
                mean_var = data_plot_mean.loc[group_i]
                std_var = data_plot_std.loc[group_i]
                x = np.arange(data_plot_mean.shape[1])
                ax.plot(mean_var, label=group_names_alt[num_group], linewidth=2, color=group_colors[num_group])
                ax.fill_between(x, mean_var - std_var, mean_var + std_var, alpha=0.2, color=group_colors[num_group])

            # --------------
            if test_type is not None:
                for group_j in group_names:
                    tract_df_name = "{}__{}_vs_{}__stats_p_vals.csv".format(tract_name, group_i, group_j)
                    tract_df_path = os.path.join(output_dir, tract_df_name)
                    inverse_tract_df_name = "{}__{}_vs_{}__stats_p_vals.csv".format(tract_name, group_j, group_i)
                    inverse_tract_df_path = os.path.join(output_dir, inverse_tract_df_name)
                    if group_i == group_j or os.path.exists(inverse_tract_df_path):
                        continue

                    print("Stats for {}".format(tract_df_name))

                    if not overwrite and os.path.exists(tract_df_path):
                        print("Skipping...already done")
                        tract_already_done = True
                        continue

                    t_test_results = compute_t_tests(column_list, tract_df, group_i, group_j, test_type=test_type,
                                                     tract_name=tract_name, correct_after_perms=correct_after_perms)

                    results_tract_df = pd.DataFrame.from_dict(t_test_results, orient='Index')
                    results_tract_df.columns = ['statistic', 'pvalue']
                    results_tract_df.insert(0, "variable", results_tract_df.index)

                    results_tract_df.to_csv(tract_df_path, index=False)

        if data_name == "afq":
            if test_type is not None and not tract_already_done:
                significant_segments_bool = results_tract_df["pvalue"] < 0.05  # 0.6
                significant_segments = results_tract_df[significant_segments_bool]
                if significant_segments.shape[0] > 0:
                    x = significant_segments["variable"].astype(int).values.tolist()
                    if "fa" in measure:
                        y = [0.2] * significant_segments.shape[0]
                    elif "tr" in measure:
                        y = [1] * significant_segments.shape[0]

                    ax.scatter(x, y, linewidth=4, marker="_", color="black")  #, label="Significant")
                    # y1 = data_plot.loc["icd"].values.tolist()
                    # y2 = data_plot.loc["noICD"].values.tolist()
                    # x = np.arange(100)
                    # ax.fill_between(x, y1, y2, where=significant_segments_bool, facecolor='green', interpolate=True)
                    # x = np.arange(data_plot.loc[group_i].shape[1])
                    # for significant in significant_segments:
                    #     ax.plot(x=significant, y=y, kind='scatter', label='significant', color='red')

            # plt.figure(figsize=(10, 2))
            # ax.legend(loc='upper right', ncol=2, fontsize=fontsize-2)
            ax.tick_params(axis='both', which='major', labelsize=fontsize-2)
            ax.tick_params(axis='both', which='minor', labelsize=fontsize-2)
            ax.xaxis.set_major_locator(plt.MaxNLocator(10))
            fig.savefig(fig_path)  #, dpi=400)
            # plt.tight_layout()
            # x1, x2, y1, y2 = plt.axis()
            # plt.axis((x1, x2, y1, y2 + 0.1))
            plt.show()


def compute_t_tests(column_list, tract_df, group_i, group_j, test_type="t", tract_name="",
                    show_tfce_plot=False, correct_after_perms=True):
    t_test_results = {}
    group1 = tract_df.where(tract_df["group"] == group_i).dropna()[column_list].values
    group2 = tract_df.where(tract_df["group"] == group_j).dropna()[column_list].values

    if test_type == "permutation":
        # https://mne.tools/stable/generated/mne.stats.permutation_cluster_test.html#mne.stats.permutation_cluster_test
        t_raw, _ = scipy.stats.ttest_ind(group1, group2, equal_var=False)
        t_raw = np.abs(t_raw)
        # ---------------------- TFCE ---------------------------
        t_tfce = cluster.tfce(t_raw, dt=.1, E=0.5, H=2, connectivity=None)

        # # ---------------------- Permutations after TFCE ---------------------------
        t_obs, p_uncorr, significant, p_clusters, cstat, ref_cstat, ref_clusters = permutation_test(group1, group2, alpha=0.05,
                                                                                      threshold=0.1,
                                                                                      n_permutations=1000,
                                                                                      equal_var=False,
                                                                                      t_vals=t_tfce,
                                                                                      cluster_stat="tmax",
                                                                                      abs=True)

        if correct_after_perms:
            significant_segments_bool, p_perm, _, _ = multipletests(pvals=p_uncorr, alpha=0.05,
                                                                    method="holm")
        else:
            p_perm = np.zeros(len(t_obs))
            for i, idx_cluster in enumerate(ref_clusters):
                p_val = 1
                if idx_cluster > 0:
                    p_val = p_clusters[idx_cluster - 1]
                p_perm[i] = p_val

        for i, column in enumerate(column_list):
            t_test_results[column] = [0, p_perm[i]]

    elif test_type == "t":
        t_real = []
        for column in column_list:
            group1 = tract_df.where(tract_df["group"] == group_i).dropna()[column]
            group2 = tract_df.where(tract_df["group"] == group_j).dropna()[column]

            t_val, p_val = scipy.stats.ttest_ind(group1, group2)
            t_test_results[column] = [t_val, p_val]
            t_val = np.abs(t_val)
            t_real.append(t_val)

    if show_tfce_plot:
        fig, ax = plt.subplots()
        plot_name = "t-values"
        ax.set_title("{} for {}".format(plot_name, tract_name))
        ax.set_ylabel(plot_name)
        ax.set_xlabel("segments along tract")
        ax.plot(t_raw, label="t original")
        ax.plot(t_tfce, label="t tfce")
        ax.legend(loc='best')
        plt.show()

    return t_test_results

def correlation_other_features(corr_features_filename, tract_var="afq", corr_features_var="all", measure="fa",
                               fontsize=16, include_hc=False):
    corr_features_data_path = os.path.join(analysis_results_dir, corr_features_filename)
    corr_features_df = pd.read_csv(corr_features_data_path)

    data_dir = get_output_dir(tract_var)
    output_dir = os.path.join(all_results_dir, data_dir)
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)

    group_dirs = [f.name for f in Path(analysis_results_dir).iterdir() if f.is_dir() and f.name.startswith("group_")]
    tract_vals_col = "tract_vals"

    all_data = {}
    for group_dir in group_dirs:
        subj_results_dir = os.path.join(analysis_results_dir, group_dir, data_dir)

        group = group_dir.split("_")[1]
        print("---- Reading {} group data".format(group))

        filenames = sorted([f.name for f in Path(subj_results_dir).iterdir() if f.suffix == ".csv" and
                            not "__stats_p_vals" in f.name and "__{}".format(measure) in f.name])
        cols_remove = ["subj_id", "group"]
        for tract_name in filenames:
            print("Reading data for {}".format(tract_name))
            tract_df_path = os.path.join(subj_results_dir, tract_name)
            tract_name = Path(tract_name).stem
            tract_df = pd.read_csv(tract_df_path, index_col=False)
            subj_ids = tract_df["subj_id"]
            corr_features_filter_df = corr_features_df.loc[:, "subj_id"].isin(subj_ids)
            corr_features_var_df = corr_features_df.loc[corr_features_filter_df]
            subj_ids_behav = corr_features_var_df["subj_id"]
            tract_filter_df = tract_df.loc[:, "subj_id"].isin(subj_ids_behav)
            tract_df = tract_df.loc[tract_filter_df]

            tract_cols = list(tract_df.columns)
            tract_cols = [elem for elem in tract_cols if elem not in cols_remove]
            x = tract_df.loc[:, tract_cols]
            x = x.mean(axis=1)
            if corr_features_var == "all":
                cols = list(corr_features_df.columns)
                cols = [elem for elem in cols if elem not in cols_remove]
                if not tract_name in all_data:
                    all_data[tract_name] = {}
                if not tract_vals_col in all_data[tract_name]:
                    all_data[tract_name][tract_vals_col] = {}
                all_data[tract_name][tract_vals_col][group] = x

                for col in cols:
                    # fig = plt.figure()
                    y = corr_features_var_df[col]
                    if not col in all_data[tract_name]:
                        all_data[tract_name][col] = {}
                    all_data[tract_name][col][group] = y

    for tract_name, vars_data in all_data.items():
        # plt.figure(figsize=(10, 10))
        x_all = vars_data[tract_vals_col]
        groups = list(x_all.keys())
        behav_vars = list(vars_data.keys())
        behav_vars.remove(tract_vals_col)

        for behav_var_name in behav_vars:
            if include_hc:
                color_groups = ["black", "#ff2600", "#73fdff"]
                group_names_alt = ["HC", "Hypersexual", "Non-hypersexual"]
            else:
                group_names_alt = ["Hypersexual", "Non-hypersexual"]
                color_groups = ["#ff2600", "#73fdff"]  # Red and turquoise
            r_groups = []
            n_groups = []
            for i, group in enumerate(groups):
                x = x_all[group]
                y = vars_data[behav_var_name][group]
                n = x.shape[0]

                try:
                    r, p = scipy.stats.pearsonr(x, y)
                    r_groups.append(r)
                    n_groups.append(n)

                    slope, intercept, r_value, p_value, std_err = scipy.stats.linregress(x, y)
                    reg_line = slope * x + intercept

                    if i >= len(color_groups):
                        raise Exception("Not enough color names to represent all groups. Add more colors names in the script")
                    plt.scatter(x, y, s=50, color=color_groups[i])  #  label=group_names_alt[i] + " (n={})".format(n)
                    plt.plot(x, reg_line, 'r', label='r={:.2f} p={:.2f}'.format(r, p), color=color_groups[i])
                except Exception as e:
                    print("Exception: ", e)

            # https://github.com/psinger/CorrelationStats/blob/master/corrstats.py
            z_fisher, p_fisher = corrstats.independent_corr(r_groups[0], r_groups[1],
                                                  n_groups[0], n_groups[1], method='fisher')

            if include_hc:
                p_fisher_label = "p_val_fisher HC-ICD corrs diff: {:.2f}".format(p_fisher)
            else:
                p_fisher_label = "p-val Fisher = {:.2f}".format(p_fisher)
            # plt.plot([], [], ' ', label=)
            plt.legend(loc='upper right', fontsize=fontsize, frameon=False)  # , ncol=2

            ax = plt.gca()
            if "DCM" in corr_features_filename:
                pass
                # ax.set_ylabel(behav_var_name, fontsize=fontsize)
            # ax.set_xlabel(tract_name, fontsize=10)
            ax.xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
            ax.tick_params(axis='both', which='major', labelsize=fontsize-2)
            ax.tick_params(axis='both', which='minor', labelsize=fontsize-2)
            ax.xaxis.set_major_locator(plt.MaxNLocator(10))
            x1, x2, y1, y2 = plt.axis()
            if "behavioral_final_med_on" in corr_features_filename:
                plt.axis((x1, x2, y1-50, y2 + 100))
            else:
                plt.axis((x1, x2, y1-0.02, y2 + 0.04))

            plt.title(p_fisher_label, fontsize=fontsize+2)
            fig_name = "corr__{}__{}".format(tract_name, behav_var_name)
            fig_path = os.path.join(all_results_corr_features_dir, fig_name)
            plt.savefig(fig_path)  #, dpi=400
            plt.show()

    return 0


# ==========================================================='

##########################################################
# ----------------- Main ------------------
##########################################################

if __name__ == "__main__":
    run_dti_estimation = False
    run_dti_normalization = False
    run_fiber_tracking = False
    run_bundle_segmentation = False
    run_get_stats = False
    run_stats_p_vals = True
    run_correlation_behavior = False

    if run_dti_estimation:
        dti_estimation_nlls(fit_method="NLLS", overwrite=False)

    if run_dti_normalization:
        tensor_normalization(preprocessing=True, force_population_mean=False, overwrite=False)

    if run_fiber_tracking:
        fiber_tracking()

    if run_bundle_segmentation:
        bundles_segmentation(run_slr=True, run_reco_bundles=True)

    if run_get_stats:
        # vba(merge_fa=True, smooth=4, num_permutations=100)
        afq_profiles(overwrite=False, show_interactive_bundle=False, measure="tr")
        bundles_shape_analysis(dsi_studio_cmd_ran=True)

    if run_stats_p_vals:
        stats_t_test("afq", only_hc_pd=False, test_type="permutation",
                     measure="tr", correct_after_perms=True, overwrite=True)  # Measure = fa OR tr
        # stats_t_test("shape")

    if run_correlation_behavior:
        corr_features_filename = "behavioral_final_med_on.csv"
        # corr_features_filename = "DCM_m12_4.csv"
        correlation_other_features(corr_features_filename, tract_var="afq", corr_features_var="all",
                                   measure="tr", include_hc=False)  # Measure=tr or fa