import os

def create_key(template, outtype=('nii.gz',), annotation_classes=None):
    if template is None or not template:
        raise ValueError('Template must be a valid format string')
    return template, outtype, annotation_classes


def infotodict(seqinfo):
    """Heuristic evaluator for determining which runs belong where
    allowed template fields - follow python string module:
    item: index within category
    subject: participant id
    seqitem: run number during scanning
    subindex: sub index within group
    """
    task_name = "david_mata_parkinson"
    t1_protocol_name = '3D_T1_MPRAGE_SAG'
    t2_bold_protocol_name = "EPI_BLOQUE"
    t2_bold_name = "bold"
    num_runs = 5
    num_run_to_zero_index = True

    t1 = create_key('sub-{subject}/{session}/anat/sub-{subject}_{session}_T1w')
    # fmap_magnitude = create_key('sub-{subject}/{session}/fmap/sub-{subject}_{session}_magnitude')
    # fmap_phase_diff = create_key('sub-{subject}/{session}/fmap/sub-{subject}_{session}_phasediff')
    dwi = create_key('sub-{subject}/{session}/dwi/sub-{subject}_{session}_run-0{item:01d}_dwi')
    fmap = create_key('sub-{subject}/{session}/fmap/sub-{subject}_{session}_dir-{dir}_epi')
    # fmap_ap_run1 = create_key('sub-{subject}/{session}/fmap/sub-{subject}_{session}_dir-AP_run-01_epi')

    t_bold_map_keys = {}
    t_bold = {}
    for i in range(num_runs):
        num_run = i if num_run_to_zero_index else i+1
        t_bold_key = create_key("sub-{subject}/{session}/func/sub-{subject}_{session}_task-" + task_name + "_run-" + str(num_run) + "_" + t2_bold_name)
        t_bold[t_bold_key] = []
        t_bold_map_keys[t2_bold_name + "_" + str(num_run)] = t_bold_key

    # t1w_high_res = create_key('sub-{subject}/{session}/anat/sub-{subject}_{session}_acq-highres_t2_bold_name + "_" + str(num_run)T1w')
    # t2w = create_key('sub-{subject}/{session}/anat/sub-{subject}_{session}_T2w')
    # flanker = create_key('sub-{subject}/{session}/func/sub-{subjalonso_m_08ect}_{session}_task-flanker_bold')
    # dwi_64dir = create_key('sub-{subject}/{session}/dwi/sub-{subject}_{session}_dwi')
    # dwi_b0s = create_key('sub-{subject}/{session}/dwi/sub-{subject}_{session}_acq-B0s_dwi')
    # rest = create_key('sub-{subject}/{session}/func/sub-{subject}_{session}_task-rest_bold')

    info = {
        t1: [],
        # fmap_magnitude: [],
        # fmap_phase_diff: [],
        dwi: [],
        fmap: [],
        # t2_bold_1: [],
        # t2_bold_2: []
        # t1w_high_res: [],
        # t2w: [],
        # flanker: [],
        # dwi_64dir: [],
        # dwi_b0s: [],
        # rest: [],
    }

    info.update(t_bold)

    for s in seqinfo:
        """
        The namedtuple `s` contains the following fields:
        * total_files_till_now
        * example_dcm_file
        * series_id
        * dcm_dir_name
        * unspecified2
        * unspecified3
        * dim1
        * dim2
        * dim3
        * dim4
        * TR
        * TE
        * protocol_name
        * is_motion_corrected
        * is_derived
        * patient_id
        * study_description
        * referring_physician_name
        * series_description
        * image_type
        """
        # info[t1] = [s.series_id]

        # print("========= RUN {} ==========".format(s.protocol_name))

        if s.protocol_name == t1_protocol_name:
            info[t1].append(s.series_id)

        # if s.series_description.startswith(t2_bold_protocol_name):
        #     num_run = int(s.series_description.split(t2_bold_protocol_name)[1])
        #     if num_run_to_zero_index:
        #         num_run -= 1
        #     print("========= BOLD {} RUN {} ==========".format(s.series_description, num_run))
        #     print("---- {} ------".format(s))
        #     t_bold_key = t_bold_map_keys[t2_bold_name + "_" + str(num_run)]
        #     info[t_bold_key].append(s.series_id)

        # if (s.dim3 == 136) and (s.dim4 == 1) and ('gre_field_mapping' in s.protocol_name):
        #     info[fmap_magnitude] = [s.series_id]
        # if (s.dim3 == 68) and (s.dim4 == 1) and ('gre_field_mapping' in s.protocol_name):
        #     info[fmap_phase_diff] = [s.series_id]

        # if s.dim3 == 70 and 'GRE_fieldmap_RS' in s.protocol_name:
        #     info[fmap_magnitude] = [s.series_id]
        # if s.dim3 == 35 and 'GRE_fieldmap_RS' in s.protocol_name:
        #     info[fmap_phase_diff] = [s.series_id]

        # if s.series_description == 'ep2d_diff_mddw_12dir_tract':
        #     info[dwi].append([s.series_id])

        if "DWI" in s.protocol_name and "dir" in s.protocol_name and \
                (not "TRACEW" in s.protocol_name) and s.dim4 != 1 :  # Spin echo fieldmap AP and PA
            info[dwi].append([s.series_id])

        if "DWI" in s.protocol_name and (not "dir" in s.protocol_name) and s.dim4 == 1:  # Spin echo fieldmap AP and PA
            dirtype = s.protocol_name.split('_')[-1]
            info[fmap].append({'item': s.series_id, 'dir': dirtype})

        print("S Image Type: ", s.image_type)
        print("S: ", s)

        # print("INFO: ", info)

    return info
