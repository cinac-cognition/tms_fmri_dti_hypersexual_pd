%fMRI behavioural onsets and volumes integration
clear; close all; clc;

for id_dataset = 1:5
    switch id_dataset
        case 1
            source_dir = 'data_path/ICD/OFF';
            struct_dir = 'data_path/ICD/Structural';
            res_dir =    'data_path/ICD/ICDOFF_SPM';
            subjs = dir([source_dir '/01*']);
        case 2
            source_dir = 'data_path/ICD/ON';
            struct_dir = 'data_path/ICD/Structural';
            res_dir =    'data_path/ICD/ICDON_SPM';
            subjs = dir([source_dir '/01*']);
        case 3
            source_dir = 'data_path/PD/ON';
            struct_dir = 'data_path/PD/Structural';
            res_dir =    'data_path/PD/PDON_SPM';
            subjs = dir([source_dir '/02*']);
        case 4
            source_dir = 'data_path/PD/OFF';
            struct_dir = 'data_path/PD/Structural';
            res_dir =    'data_path/PD/PDOFF_SPM';
            subjs = dir([source_dir '/02*']);
        case 5
            source_dir = 'data_path/HC';
            struct_dir = 'data_path/HC/Structural';
            res_dir =    'data_path/HC';
            subjs = dir([source_dir '/03*']);
    end
    
    if ~exist(res_dir,'dir'); mkdir(res_dir); end
    fslpreline = 'export FSLDIR=/usr/local/fsl ; source $FSLDIR/etc/fslconf/fsl.sh ; $FSLDIR/bin/';
    
    prefix_dcm2nii = 'vol_';
    
    % Fieldmap params (same for all subjects)
    param.esp = 0.53;
    param.asset = 0.5;
    
    PATHS.spm12_path = '/path/spm12/';
    rmpath(genpath(PATHS.spm12_path));
    addpath(genpath(PATHS.spm12_path));
    
    for id_s = 1:length(subjs)
        
        do_fmap = 1;
        subjname = subjs(id_s).name;
        subject_dir = fullfile(res_dir,subjname);
        
        % Preproc T1
        t1folderin = fullfile(struct_dir,subjname);
        t1dcmfolder = dir([t1folderin '/3D*']);
        T1dir = fullfile(t1folderin,'T1');
        if ~exist(T1dir,'dir'); mkdir(T1dir); end
        T1file = [t1folderin '/' t1dcmfolder(1).name '/*.nii'];
        pathT1 = fullfile(T1dir,'T1.nii');
        fmrifolderin = fullfile(res_dir,subjname,'fMRI_all');
        fmrifolderin_spmstats = fullfile(res_dir,subjname,'fMRI_all_Stat');
        if ~exist(fmrifolderin_spmstats,'dir'); mkdir(fmrifolderin_spmstats); end
        
        %------------------------------------------------------------------
        % Run new_segment + DARTEL
        if ~exist(fullfile(T1dir,'rc3T1.nii'),'file')
            clear matlabbatch; spm_jobman('initcfg');
            matlabbatch{1}.spm.spatial.preproc.channel.vols = {[T1dir '/T1.nii,1']};
            matlabbatch{1}.spm.spatial.preproc.channel.biasreg = 0.001;
            matlabbatch{1}.spm.spatial.preproc.channel.biasfwhm = 60;
            matlabbatch{1}.spm.spatial.preproc.channel.write = [0 0];
            matlabbatch{1}.spm.spatial.preproc.tissue(1).tpm = {[PATHS.spm12_path '/tpm/TPM.nii,1']};
            matlabbatch{1}.spm.spatial.preproc.tissue(1).ngaus = 1;
            matlabbatch{1}.spm.spatial.preproc.tissue(1).native = [1 1];
            matlabbatch{1}.spm.spatial.preproc.tissue(1).warped = [0 0];
            matlabbatch{1}.spm.spatial.preproc.tissue(2).tpm = {[PATHS.spm12_path '/tpm/TPM.nii,2']};
            matlabbatch{1}.spm.spatial.preproc.tissue(2).ngaus = 1;
            matlabbatch{1}.spm.spatial.preproc.tissue(2).native = [1 1];
            matlabbatch{1}.spm.spatial.preproc.tissue(2).warped = [0 1];
            matlabbatch{1}.spm.spatial.preproc.tissue(3).tpm = {[PATHS.spm12_path '/tpm/TPM.nii,3']};
            matlabbatch{1}.spm.spatial.preproc.tissue(3).ngaus = 2;
            matlabbatch{1}.spm.spatial.preproc.tissue(3).native = [1 1];
            matlabbatch{1}.spm.spatial.preproc.tissue(3).warped = [0 0];
            matlabbatch{1}.spm.spatial.preproc.tissue(4).tpm = {[PATHS.spm12_path '/tpm/TPM.nii,4']};
            matlabbatch{1}.spm.spatial.preproc.tissue(4).ngaus = 3;
            matlabbatch{1}.spm.spatial.preproc.tissue(4).native = [0 0];
            matlabbatch{1}.spm.spatial.preproc.tissue(4).warped = [0 0];
            matlabbatch{1}.spm.spatial.preproc.tissue(5).tpm = {[PATHS.spm12_path '/tpm/TPM.nii,5']};
            matlabbatch{1}.spm.spatial.preproc.tissue(5).ngaus = 4;
            matlabbatch{1}.spm.spatial.preproc.tissue(5).native = [0 0];
            matlabbatch{1}.spm.spatial.preproc.tissue(5).warped = [0 0];
            matlabbatch{1}.spm.spatial.preproc.tissue(6).tpm = {[PATHS.spm12_path '/tpm/TPM.nii,6']};
            matlabbatch{1}.spm.spatial.preproc.tissue(6).ngaus = 2;
            matlabbatch{1}.spm.spatial.preproc.tissue(6).native = [0 0];
            matlabbatch{1}.spm.spatial.preproc.tissue(6).warped = [0 0];
            matlabbatch{1}.spm.spatial.preproc.warp.mrf = 1;
            matlabbatch{1}.spm.spatial.preproc.warp.cleanup = 1;
            matlabbatch{1}.spm.spatial.preproc.warp.reg = [0 0.001 0.5 0.05 0.2];
            matlabbatch{1}.spm.spatial.preproc.warp.affreg = 'mni';
            matlabbatch{1}.spm.spatial.preproc.warp.fwhm = 0;
            matlabbatch{1}.spm.spatial.preproc.warp.samp = 3;
            matlabbatch{1}.spm.spatial.preproc.warp.write = [1 1];
            spm('defaults', 'PET');
            spm_jobman('run', matlabbatch);
        end        
        
        %------------------------------------------------------------------
        % Fieldmap preparation
        clc;
        fprintf('........................\n');
        fprintf('Fieldmap preparation.\n');
        
        subjdir = fullfile(source_dir,subjname);
        fmap_folders = dir([subjdir '/GRE_*']);
        fmap_dir = fullfile(subject_dir,'fmap');
        if ~exist(fmap_dir,'dir'); mkdir(fmap_dir); end
        
        if ~isempty(fmap_folders)
            % Convert to nifti magnitud image and fieldmap
            fmap_dir_dcm = [subjdir '/' fmap_folders(end).name];
            runline = ['dcm2nii -o ' fmap_dir ' ' fmap_dir_dcm];
            [status,res] = system(runline);
            system(['gzip ' fmap_dir '/GRE_*nii']);
            nfmapfiles = dir([fmap_dir '/GRE_*nii.gz']);
            if length(nfmapfiles) > 1
                delete([fmap_dir '/' nfmapfiles(1).name]);
            end
            [status,res] = system(['mv ' fmap_dir '/GRE_*nii.gz ' fmap_dir '/fmap.nii.gz']);
            
            mag_dir_dcm = [subjdir '/' fmap_folders(1).name];
            runline = ['dcm2nii -o ' fmap_dir ' ' mag_dir_dcm];
            [status,res] = system(runline);
            system(['gzip ' fmap_dir '/GRE_*nii']);
            nfmapfiles = dir([fmap_dir '/GRE_*nii.gz']);
            if length(nfmapfiles) > 1
                delete([fmap_dir '/' nfmapfiles(1).name]);
            end
            [status,res] = system(['mv ' fmap_dir '/GRE_*nii.gz ' fmap_dir '/mag.nii.gz']);
            
            % Bet Magnitude Image
            fileIn = fullfile(fmap_dir,'mag.nii.gz');
            fileOut = fullfile(fmap_dir,'mag_brain.nii.gz');
            sentence = sprintf('%sbet %s %s -f 0.3 -m',fslpreline,fileIn,fileOut);
            [status,result] = system(sentence);
            % Brain mask fieldmap
            fileIn = fullfile(fmap_dir,'fmap.nii.gz');
            fileMas = fullfile(fmap_dir,'mag_brain_mask.nii.gz');
            fileOut = fullfile(fmap_dir,'fmap.nii.gz');
            sentence = sprintf('%sfslmaths %s -mas %s %s',fslpreline,fileIn,fileMas,fileOut);
            [status,result] = system(sentence);
            
            %----------------------------------------------------------------------
            % Prepare fieldmap
            fileIn = fullfile(fmap_dir,'fmap.nii.gz');
            fileMagRes = fullfile(fmap_dir,'mag_brain.nii.gz');
            fileOut = fullfile(fmap_dir,'fmap_rads.nii.gz');
            sentence = sprintf('%sfsl_prepare_fieldmap SIEMENS %s %s %s %f',fslpreline,fileIn,fileMagRes,fileOut,2.46);
            [status,result] = system(sentence);
            
            % Reslice fmap into fMRI dimensions
            fileRef = fullfile(fmrifolderin,'vol_T1_0001.nii');
            fileIn = fullfile(fmap_dir,'fmap_rads.nii.gz');
            fileInRes = fullfile(fmap_dir,'fmap_rads.nii.gz');
            sentence = sprintf('WarpImageMultiTransform 3 %s %s -R %s --tightest-bounding-box --reslice-by-header',fileIn,fileInRes,fileRef);
            [status,result] = system(sentence);
            fileIn = fullfile(fmap_dir,'mag_brain.nii.gz');
            fileInRes = fullfile(fmap_dir,'mag_brain.nii.gz');
            sentence = sprintf('WarpImageMultiTransform 3 %s %s -R %s --tightest-bounding-box --reslice-by-header',fileIn,fileInRes,fileRef);
            [status,result] = system(sentence);
            
            fileIn = fullfile(fmap_dir,'mag_brain_mask.nii.gz');
            fileInRes = fullfile(fmap_dir,'mag_brain_mask.nii.gz');
            sentence = sprintf('WarpImageMultiTransform 3 %s %s -R %s --tightest-bounding-box --reslice-by-header',fileIn,fileInRes,fileRef);
            [status,result] = system(sentence);
            
            % Delete residual files
            file1 = fullfile(fmap_dir,'mag.nii.gz');
            delete(file1);
        else
            do_fmap = 0;
        end
        
        %--------------------------------------------------------------------------
        % Prepare files and folders in fMRI data
        fprintf('........................\n');
        fprintf('Preparing paths, and running slice timing and realignment.\n');
        clear matlabbatch; spm_jobman('initcfg');
        epivol = dir([fmrifolderin '/vol_T*']);
        clear paths_vols;
        for id_epi = 1:length(epivol)
            paths_vols{id_epi,1} = [fmrifolderin '/' epivol(id_epi).name ',1'];
        end
        param.NSlices = 35;
        param.TR = 2000;
        
        matlabbatch{1}.spm.temporal.st.scans = {paths_vols}';
        matlabbatch{1}.spm.temporal.st.nslices = param.NSlices;
        matlabbatch{1}.spm.temporal.st.tr = param.TR/1000;
        matlabbatch{1}.spm.temporal.st.ta = (param.TR/1000) - ((param.TR/1000)/(param.NSlices));
        matlabbatch{1}.spm.temporal.st.so = [1:2:param.NSlices 2:2:param.NSlices];
        matlabbatch{1}.spm.temporal.st.refslice = 1;
        matlabbatch{1}.spm.temporal.st.prefix = 'a';
        spm('defaults', 'FMRI');
        spm_jobman('run', matlabbatch);
        clear paths_vols paths_vols1 paths_vols2 paths_vols3 paths_vols4 paths_vols5;
        clear matlabbatch;
        spm_jobman('initcfg');
        volumes1 = dir([fmrifolderin '/avol_T1*']);
        for mv = 1:length(volumes1)
            paths_vols1{mv,1} = [fmrifolderin '/' volumes1(mv).name ',1'];
        end
        volumes2 = dir([fmrifolderin '/avol_T2*']);
        for mv = 1:length(volumes2)
            paths_vols2{mv,1} = [fmrifolderin '/' volumes2(mv).name ',1'];
        end
        volumes3 = dir([fmrifolderin '/avol_T3*']);
        for mv = 1:length(volumes3)
            paths_vols3{mv,1} = [fmrifolderin '/' volumes3(mv).name ',1'];
        end
        volumes4 = dir([fmrifolderin '/avol_T4*']);
        for mv = 1:length(volumes4)
            paths_vols4{mv,1} = [fmrifolderin '/' volumes4(mv).name ',1'];
        end
        volumes5 = dir([fmrifolderin '/avol_T5*']);
        if ~isempty(volumes5)
            for mv = 1:length(volumes5)
                paths_vols5{mv,1} = [fmrifolderin '/' volumes5(mv).name ',1'];
            end
            matlabbatch{1}.spm.spatial.realign.estwrite.data = {paths_vols1
                paths_vols2
                paths_vols3
                paths_vols4
                paths_vols5}';
        else
            matlabbatch{1}.spm.spatial.realign.estwrite.data = {paths_vols1
                paths_vols2
                paths_vols3
                paths_vols4}';
        end
        matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.quality = 0.9;
        matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.sep = 4;
        matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.fwhm = 5;
        matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.rtm = 1;
        matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.interp = 2;
        matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.wrap = [0 0 0];
        matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.weight = '';
        matlabbatch{1}.spm.spatial.realign.estwrite.roptions.which = [2 1];
        matlabbatch{1}.spm.spatial.realign.estwrite.roptions.interp = 4;
        matlabbatch{1}.spm.spatial.realign.estwrite.roptions.wrap = [0 0 0];
        matlabbatch{1}.spm.spatial.realign.estwrite.roptions.mask = 1;
        matlabbatch{1}.spm.spatial.realign.estwrite.roptions.prefix = 'r';
        spm('defaults', 'FMRI');
        spm_jobman('run', matlabbatch);
        clear paths_vols1 paths_vols2 paths_vols3 paths_vols4 paths_vols5;
        
        %------------------------------------------------------------------
        % Fieldmap correction for EPI
        fprintf('........................\n');
        fprintf('Running fieldmap correction.\n');
        subjdir = fullfile(source_dir,subjname);
        fmap_dir = fullfile(subject_dir,'fmap');
        
        if do_fmap
            
            %------------------------------------------------------------------
            % Unwarp EPI data
            fileOut = fullfile(fmap_dir,'merged4Depi.nii.gz');
            sentence = sprintf('%sfslmerge -t %s %s/ravol*',fslpreline,fileOut,fmrifolderin);
            [status,result] = system(sentence);
            
            fileIn = fullfile(fmap_dir,'merged4Depi.nii.gz');
            filefmap = fullfile(fmap_dir,'fmap_rads.nii.gz');
            fileout = fullfile(fmap_dir,'merged4Depi_dewarped.nii.gz');
            fileMask = fullfile(fmap_dir,'mag_brain_mask.nii.gz');
            
            runline = sprintf('%sfugue -i %s --dwell=%f --loadfmap=%s -u %s --unwarpdir=y- --mask=%s -s 0.5',fslpreline,fileIn,(param.esp*param.asset)/1000,filefmap,fileout,fileMask);
            [status,result] = system(runline);
            delete(fileIn);
            pause(30);
            %------------------------------------------------------------------
            % Extract mean dewarped volume as a reference
            fileIn = fullfile(fmap_dir,'merged4Depi_dewarped.nii.gz');
            fileOut = fullfile(fmap_dir,'exf_brain.nii.gz');
            sentence = sprintf('%sfslmaths %s -Tmean %s',fslpreline,fileIn,fileOut);
            [status,result] = system(sentence);
            
            %------------------------------------------------------------------
            % Get Mean fMRI, brain extract, then convert to nifti fieldmap,
            % brain extract, coregister to fMRI and reslice, prepare fieldmap
            % with fsl and apply fieldmap correction
            fileIn = fullfile(fmap_dir,'exf_brain.nii.gz');
            fileOut = fullfile(fmrifolderin,'EPI_mean.nii.gz');
            system(['mv ' fileIn ' ' fileOut]);
            system(['gunzip ' fileOut]);
            fileIn = fullfile(fmap_dir,'merged4Depi_dewarped.nii.gz');
            fileOut = fullfile(fmrifolderin,'dravol');
            sentence = sprintf('%sfslsplit %s %s/dravol -t',fslpreline,fileIn,fmrifolderin);
            [status,result] = system(sentence);
            system(['gunzip ' fmrifolderin '/dravol*']);
            delete(fileIn);
        else
            
            fileIn = fullfile(fmap_dir,'merged4Depi.nii.gz');
            sentence = sprintf('%sfslmerge -t %s %s/ravol*',fslpreline,fileIn,fmrifolderin);
            [status,result] = system(sentence);
            
            fileOut = fullfile(fmap_dir,'exf_brain.nii.gz');
            sentence = sprintf('%sfslmaths %s -Tmean %s',fslpreline,fileIn,fileOut);
            [status,result] = system(sentence);
            
            fileIn = fullfile(fmap_dir,'exf_brain.nii.gz');
            fileOut = fullfile(fmrifolderin,'EPI_mean.nii.gz');
            system(['mv ' fileIn ' ' fileOut]);
            system(['gunzip ' fileOut]);
            
            fileIn = fullfile(fmap_dir,'merged4Depi.nii.gz');
            fileOut = fullfile(fmrifolderin,'dravol');
            sentence = sprintf('%sfslsplit %s %s/dravol -t',fslpreline,fileIn,fmrifolderin);
            [status,result] = system(sentence);
            system(['gunzip ' fmrifolderin '/dravol*']);
            delete(fileIn);
            
        end
        %------------------------------------------------------------------
        % Estimate coregister to T1 high resolution image, not reslice
        Ref_file = fullfile(T1dir,'T1.nii');
        Source_file = fullfile(fmrifolderin,'EPI_mean.nii');
        Def_file = fullfile(T1dir,'y_T1.nii');
        
        volumes = dir([fmrifolderin '/dravol*']);
        clear paths_vols;
        for mv = 1:length(volumes)
            paths_vols{mv,1} = [fmrifolderin '/' volumes(mv).name ',1'];
        end
        clear matlabbatch; spm_jobman('initcfg');
        matlabbatch{1}.spm.spatial.coreg.estimate.ref = {[Ref_file ',1']};
        matlabbatch{1}.spm.spatial.coreg.estimate.source = {[Source_file ',1']};
        matlabbatch{1}.spm.spatial.coreg.estimate.other = paths_vols;
        matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.cost_fun = 'nmi';
        matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.sep = [4 2];
        matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
        matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.fwhm = [7 7];
        spm('defaults', 'FMRI');
        spm_jobman('run', matlabbatch);
        
        clear matlabbatch; spm_jobman('initcfg');
        matlabbatch{1}.spm.spatial.normalise.write.subj.def = {Def_file};
        matlabbatch{1}.spm.spatial.normalise.write.subj.resample = paths_vols;
        matlabbatch{1}.spm.spatial.normalise.write.woptions.bb = [NaN NaN NaN
            NaN NaN NaN];
        matlabbatch{1}.spm.spatial.normalise.write.woptions.vox = [2 2 2];
        matlabbatch{1}.spm.spatial.normalise.write.woptions.interp = 4;
        matlabbatch{2}.spm.spatial.smooth.data(1) = cfg_dep('Normalise: Write: Normalised Images (Subj 1)', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('()',{1}, '.','files'));
        matlabbatch{2}.spm.spatial.smooth.fwhm = [8 8 8];
        matlabbatch{2}.spm.spatial.smooth.dtype = 0;
        matlabbatch{2}.spm.spatial.smooth.im = 0;
        matlabbatch{2}.spm.spatial.smooth.prefix = 's';
        spm('defaults', 'FMRI');
        spm_jobman('run', matlabbatch);
        
        fileIn = fullfile(fmrifolderin,'ravol*.nii'); system(['rm ' fileIn]);
        fileIn = fullfile(fmrifolderin,'dravol*.nii'); system(['rm ' fileIn]);
        fileIn = fullfile(fmrifolderin,'avol*.nii'); system(['rm ' fileIn]);
        fileIn = fullfile(fmrifolderin,'wdravol*.nii'); system(['rm ' fileIn]);
        clear paths_vols; volumes = dir([fmrifolderin '/swdravol*']);
        
        for mv = 1:length(volumes)
            paths_vols{mv,1} = [fmrifolderin '/' volumes(mv).name ',1'];
        end
        realign_files = dir([fmrifolderin '/rp*txt']);
        
        % Derivatives and quadratic forms
        for mv = 1:length(realign_files)
            tmp = importdata([fmrifolderin '/' realign_files(mv).name]);
            if mv == 1
                tmp2 = [zeros(1,6);diff(tmp)];
                R = [tmp tmp2 tmp.^2];
            else
                tmp2 = [zeros(1,6);diff(tmp)];
                R(end+1:end+size(tmp,1),:) = [tmp tmp2 tmp.^2];
            end
        end
        
        d_circ = 2*pi*50 ; % Head is a sphere of radius equal to 50mm
        Angle2Dist = d_circ * R(:,4:6) ./ pi;
        diffR = [zeros(1,6);diff([R(:,1:3) Angle2Dist])];
        FD = sum(abs(diffR),2);       
        save([fmrifolderin '/motionres.mat'],'R');        
        
        load(fullfile(fmrifolderin,'Onsets.mat'));
        clear matlabbatch; spm_jobman('initcfg');
        matlabbatch{1}.spm.stats.fmri_spec.dir = {fmrifolderin_spmstats};
        matlabbatch{1}.spm.stats.fmri_spec.timing.units = 'secs';
        matlabbatch{1}.spm.stats.fmri_spec.timing.RT = 2;
        matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t = 16;
        matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t0 = 1;
        matlabbatch{1}.spm.stats.fmri_spec.sess.scans = paths_vols;
        matlabbatch{1}.spm.stats.fmri_spec.sess.cond(1).name = 'FixOnset';
        matlabbatch{1}.spm.stats.fmri_spec.sess.cond(1).onset = FixOnset;
        matlabbatch{1}.spm.stats.fmri_spec.sess.cond(1).duration = 0;
        matlabbatch{1}.spm.stats.fmri_spec.sess.cond(1).tmod = 0;
        matlabbatch{1}.spm.stats.fmri_spec.sess.cond(1).pmod = struct('name', {}, 'param', {}, 'poly', {});
        
        matlabbatch{1}.spm.stats.fmri_spec.sess.cond(2).name = 'StopOnsetOk_ero';
        matlabbatch{1}.spm.stats.fmri_spec.sess.cond(2).onset = StopOnsetOk_ero;
        matlabbatch{1}.spm.stats.fmri_spec.sess.cond(2).duration = 0;
        matlabbatch{1}.spm.stats.fmri_spec.sess.cond(2).tmod = 0;
        matlabbatch{1}.spm.stats.fmri_spec.sess.cond(2).pmod = struct('name', {}, 'param', {}, 'poly', {});
        
        matlabbatch{1}.spm.stats.fmri_spec.sess.cond(3).name = 'StopOnsetOk_noero';
        matlabbatch{1}.spm.stats.fmri_spec.sess.cond(3).onset = StopOnsetOk_noero;
        matlabbatch{1}.spm.stats.fmri_spec.sess.cond(3).duration = 0;
        matlabbatch{1}.spm.stats.fmri_spec.sess.cond(3).tmod = 0;
        matlabbatch{1}.spm.stats.fmri_spec.sess.cond(3).pmod = struct('name', {}, 'param', {}, 'poly', {});
        
        if isempty(StopOnsetFail_ero); StopOnsetFail_ero = [Inf    1]; end
        matlabbatch{1}.spm.stats.fmri_spec.sess.cond(4).name = 'StopOnsetFail_ero';
        matlabbatch{1}.spm.stats.fmri_spec.sess.cond(4).onset = StopOnsetFail_ero;
        matlabbatch{1}.spm.stats.fmri_spec.sess.cond(4).duration = 0;
        matlabbatch{1}.spm.stats.fmri_spec.sess.cond(4).tmod = 0;
        matlabbatch{1}.spm.stats.fmri_spec.sess.cond(4).pmod = struct('name', {}, 'param', {}, 'poly', {});
        
        if isempty(StopOnsetFail_noero); StopOnsetFail_noero = [Inf    1]; end
        matlabbatch{1}.spm.stats.fmri_spec.sess.cond(5).name = 'StopOnsetFail_noero';
        matlabbatch{1}.spm.stats.fmri_spec.sess.cond(5).onset = StopOnsetFail_noero;
        matlabbatch{1}.spm.stats.fmri_spec.sess.cond(5).duration = 0;
        matlabbatch{1}.spm.stats.fmri_spec.sess.cond(5).tmod = 0;
        matlabbatch{1}.spm.stats.fmri_spec.sess.cond(5).pmod = struct('name', {}, 'param', {}, 'poly', {});
        
        matlabbatch{1}.spm.stats.fmri_spec.sess.cond(6).name = 'PicOnset_ero';
        matlabbatch{1}.spm.stats.fmri_spec.sess.cond(6).onset = PicOnset_ero;
        matlabbatch{1}.spm.stats.fmri_spec.sess.cond(6).duration = 0;
        matlabbatch{1}.spm.stats.fmri_spec.sess.cond(6).tmod = 0;
        matlabbatch{1}.spm.stats.fmri_spec.sess.cond(6).pmod = struct('name', {}, 'param', {}, 'poly', {});
        
        matlabbatch{1}.spm.stats.fmri_spec.sess.cond(7).name = 'PicOnset_noero';
        matlabbatch{1}.spm.stats.fmri_spec.sess.cond(7).onset = PicOnset_noero;
        matlabbatch{1}.spm.stats.fmri_spec.sess.cond(7).duration = 0;
        matlabbatch{1}.spm.stats.fmri_spec.sess.cond(7).tmod = 0;
        matlabbatch{1}.spm.stats.fmri_spec.sess.cond(7).pmod = struct('name', {}, 'param', {}, 'poly', {});
        
        matlabbatch{1}.spm.stats.fmri_spec.sess.cond(8).name = 'GoOnset_ero';
        matlabbatch{1}.spm.stats.fmri_spec.sess.cond(8).onset = GoOnset_ero;
        matlabbatch{1}.spm.stats.fmri_spec.sess.cond(8).duration = 0;
        matlabbatch{1}.spm.stats.fmri_spec.sess.cond(8).tmod = 0;
        matlabbatch{1}.spm.stats.fmri_spec.sess.cond(8).pmod = struct('name', {}, 'param', {}, 'poly', {});
        
        matlabbatch{1}.spm.stats.fmri_spec.sess.cond(9).name = 'GoOnset_noero';
        matlabbatch{1}.spm.stats.fmri_spec.sess.cond(9).onset = GoOnset_noero;
        matlabbatch{1}.spm.stats.fmri_spec.sess.cond(9).duration = 0;
        matlabbatch{1}.spm.stats.fmri_spec.sess.cond(9).tmod = 0;
        matlabbatch{1}.spm.stats.fmri_spec.sess.cond(9).pmod = struct('name', {}, 'param', {}, 'poly', {});
        
        matlabbatch{1}.spm.stats.fmri_spec.sess.multi = {''};
        matlabbatch{1}.spm.stats.fmri_spec.sess.regress = struct('name', {}, 'val', {});
        
        matlabbatch{1}.spm.stats.fmri_spec.sess.multi_reg = {[fmrifolderin '/motionres.mat']};
        matlabbatch{1}.spm.stats.fmri_spec.sess.hpf = 128;
        matlabbatch{1}.spm.stats.fmri_spec.fact = struct('name', {}, 'levels', {});
        matlabbatch{1}.spm.stats.fmri_spec.bases.hrf.derivs = [1 0];
        matlabbatch{1}.spm.stats.fmri_spec.volt = 1;
        matlabbatch{1}.spm.stats.fmri_spec.global = 'None';
        matlabbatch{1}.spm.stats.fmri_spec.mthresh = 0.8;
        matlabbatch{1}.spm.stats.fmri_spec.mask = {''};
        matlabbatch{1}.spm.stats.fmri_spec.cvi = 'AR(1)';
        
        matlabbatch{2}.spm.stats.fmri_est.spmmat(1) = cfg_dep('fMRI model specification: SPM.mat File', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
        matlabbatch{2}.spm.stats.fmri_est.write_residuals = 0;
        matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;
        matlabbatch{3}.spm.stats.con.spmmat(1) = cfg_dep('Model estimation: SPM.mat File', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
        matlabbatch{3}.spm.stats.con.consess{1}.tcon.name = 'StopOK_ero';
        matlabbatch{3}.spm.stats.con.consess{1}.tcon.weights = [0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];  
        matlabbatch{3}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
        matlabbatch{3}.spm.stats.con.consess{2}.tcon.name = 'StopOK_Nero';
        matlabbatch{3}.spm.stats.con.consess{2}.tcon.weights = [0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0];  
        matlabbatch{3}.spm.stats.con.consess{2}.tcon.sessrep = 'none';
        matlabbatch{3}.spm.stats.con.consess{3}.tcon.name = 'StopOK_ero > StopOK_Nero';
        matlabbatch{3}.spm.stats.con.consess{3}.tcon.weights = [0 0 1 0 -1 0 0 0 0 0 0 0 0 0 0 0 0 0];  
        matlabbatch{3}.spm.stats.con.consess{3}.tcon.sessrep = 'none';
        matlabbatch{3}.spm.stats.con.consess{4}.tcon.name = 'StopOK_Nero > StopOK_ero';
        matlabbatch{3}.spm.stats.con.consess{4}.tcon.weights = [0 0 -1 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0];  
        matlabbatch{3}.spm.stats.con.consess{4}.tcon.sessrep = 'none';
        matlabbatch{3}.spm.stats.con.consess{5}.tcon.name = 'StopOK_ero > Null';
        matlabbatch{3}.spm.stats.con.consess{5}.tcon.weights = [-1 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];  
        matlabbatch{3}.spm.stats.con.consess{5}.tcon.sessrep = 'none';
        matlabbatch{3}.spm.stats.con.consess{6}.tcon.name = 'StopOK_Nero > Null';
        matlabbatch{3}.spm.stats.con.consess{6}.tcon.weights = [-1 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0];  
        matlabbatch{3}.spm.stats.con.consess{6}.tcon.sessrep = 'none';
        matlabbatch{3}.spm.stats.con.consess{7}.tcon.name = 'StopOK_ero > GoEro';
        matlabbatch{3}.spm.stats.con.consess{7}.tcon.weights = [0 0 1 0 0 0 0 0 0 0 0 0 0 0 -1 0 0 0];  
        matlabbatch{3}.spm.stats.con.consess{7}.tcon.sessrep = 'none';
        matlabbatch{3}.spm.stats.con.consess{8}.tcon.name = 'StopOK_Nero > GoNEro';
        matlabbatch{3}.spm.stats.con.consess{8}.tcon.weights = [0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 -1 0];  
        matlabbatch{3}.spm.stats.con.consess{8}.tcon.sessrep = 'none';
        matlabbatch{3}.spm.stats.con.consess{9}.tcon.name = 'StopFail_ero';
        matlabbatch{3}.spm.stats.con.consess{9}.tcon.weights = [0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0];
        matlabbatch{3}.spm.stats.con.consess{9}.tcon.sessrep = 'none';
        matlabbatch{3}.spm.stats.con.consess{10}.tcon.name = 'StopFail_Nero';
        matlabbatch{3}.spm.stats.con.consess{10}.tcon.weights = [0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0];
        matlabbatch{3}.spm.stats.con.consess{10}.tcon.sessrep = 'none';
        matlabbatch{3}.spm.stats.con.consess{11}.tcon.name = 'StopFail_ero > StopFail_Nero';
        matlabbatch{3}.spm.stats.con.consess{11}.tcon.weights = [0 0 0 0 0 0 1 0 -1 0 0 0 0 0 0 0 0 0];
        matlabbatch{3}.spm.stats.con.consess{11}.tcon.sessrep = 'none';
        matlabbatch{3}.spm.stats.con.consess{12}.tcon.name = 'StopFail_Nero > StopFail_ero';
        matlabbatch{3}.spm.stats.con.consess{12}.tcon.weights = [0 0 0 0 0 0 -1 0 1 0 0 0 0 0 0 0 0 0];
        matlabbatch{3}.spm.stats.con.consess{12}.tcon.sessrep = 'none';
        matlabbatch{3}.spm.stats.con.consess{13}.tcon.name = 'StopFail_ero > Null';
        matlabbatch{3}.spm.stats.con.consess{13}.tcon.weights = [-1 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0];
        matlabbatch{3}.spm.stats.con.consess{13}.tcon.sessrep = 'none';
        matlabbatch{3}.spm.stats.con.consess{14}.tcon.name = 'StopFail_Nero > Null';
        matlabbatch{3}.spm.stats.con.consess{14}.tcon.weights = [-1 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0];
        matlabbatch{3}.spm.stats.con.consess{14}.tcon.sessrep = 'none';
        matlabbatch{3}.spm.stats.con.consess{15}.tcon.name = 'StopFail_ero > GoEro';
        matlabbatch{3}.spm.stats.con.consess{15}.tcon.weights = [0 0 0 0 0 0 1 0 0 0 0 0 0 0 -1 0 0 0];
        matlabbatch{3}.spm.stats.con.consess{15}.tcon.sessrep = 'none';
        matlabbatch{3}.spm.stats.con.consess{16}.tcon.name = 'StopFail_Nero > GoNero';
        matlabbatch{3}.spm.stats.con.consess{16}.tcon.weights = [0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 -1 0];
        matlabbatch{3}.spm.stats.con.consess{16}.tcon.sessrep = 'none';
        matlabbatch{3}.spm.stats.con.consess{17}.tcon.name = 'Image_ero';
        matlabbatch{3}.spm.stats.con.consess{17}.tcon.weights = [0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0];
        matlabbatch{3}.spm.stats.con.consess{17}.tcon.sessrep = 'none';
        matlabbatch{3}.spm.stats.con.consess{18}.tcon.name = 'Image_Nero';
        matlabbatch{3}.spm.stats.con.consess{18}.tcon.weights = [0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0];  
        matlabbatch{3}.spm.stats.con.consess{18}.tcon.sessrep = 'none';
        matlabbatch{3}.spm.stats.con.consess{19}.tcon.name = 'Image_Ero > Nero';
        matlabbatch{3}.spm.stats.con.consess{19}.tcon.weights = [0 0 0 0 0 0 0 0 0 0 1 0 -1 0 0 0 0 0];  
        matlabbatch{3}.spm.stats.con.consess{19}.tcon.sessrep = 'none';
        matlabbatch{3}.spm.stats.con.consess{20}.tcon.name = 'Image_Nero > ero';
        matlabbatch{3}.spm.stats.con.consess{20}.tcon.weights = [0 0 0 0 0 0 0 0 0 0 -1 0 1 0 0 0 0 0];  
        matlabbatch{3}.spm.stats.con.consess{20}.tcon.sessrep = 'none';
        matlabbatch{3}.spm.stats.con.consess{21}.tcon.name = 'Image_Ero > Null';
        matlabbatch{3}.spm.stats.con.consess{21}.tcon.weights = [-1 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0];  
        matlabbatch{3}.spm.stats.con.consess{21}.tcon.sessrep = 'none';
        matlabbatch{3}.spm.stats.con.consess{22}.tcon.name = 'Image_Nero > Null';
        matlabbatch{3}.spm.stats.con.consess{22}.tcon.weights = [-1 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0];  
        matlabbatch{3}.spm.stats.con.consess{22}.tcon.sessrep = 'none';
        matlabbatch{3}.spm.stats.con.consess{23}.tcon.name = 'GoEro';
        matlabbatch{3}.spm.stats.con.consess{23}.tcon.weights = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0];  
        matlabbatch{3}.spm.stats.con.consess{23}.tcon.sessrep = 'none';
        matlabbatch{3}.spm.stats.con.consess{24}.tcon.name = 'GoNEro';
        matlabbatch{3}.spm.stats.con.consess{24}.tcon.weights = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0];  
        matlabbatch{3}.spm.stats.con.consess{24}.tcon.sessrep = 'none';
        matlabbatch{3}.spm.stats.con.consess{25}.tcon.name = 'GoEro > GoNero';
        matlabbatch{3}.spm.stats.con.consess{25}.tcon.weights = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 -1 0];  
        matlabbatch{3}.spm.stats.con.consess{25}.tcon.sessrep = 'none';
        matlabbatch{3}.spm.stats.con.consess{26}.tcon.name = 'GoNEro > Goero';
        matlabbatch{3}.spm.stats.con.consess{26}.tcon.weights = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 -1 0 1 0];  
        matlabbatch{3}.spm.stats.con.consess{26}.tcon.sessrep = 'none';
        matlabbatch{3}.spm.stats.con.consess{27}.tcon.name = 'StopOnsetOk_both';
        matlabbatch{3}.spm.stats.con.consess{27}.tcon.weights = [0 0 1 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0];  
        matlabbatch{3}.spm.stats.con.consess{27}.tcon.sessrep = 'none';
        matlabbatch{3}.spm.stats.con.consess{28}.tcon.name = 'StopOnsetFail_both';
        matlabbatch{3}.spm.stats.con.consess{28}.tcon.weights = [0 0 0 0 0 0 1 0 1 0 0 0 0 0 0 0 0 0];
        matlabbatch{3}.spm.stats.con.consess{28}.tcon.sessrep = 'none';
        matlabbatch{3}.spm.stats.con.consess{29}.tcon.name = 'GoOnset_both';
        matlabbatch{3}.spm.stats.con.consess{29}.tcon.weights = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 1 0];
        matlabbatch{3}.spm.stats.con.consess{29}.tcon.sessrep = 'none';
        matlabbatch{3}.spm.stats.con.consess{30}.tcon.name = 'ImageOnset_both';
        matlabbatch{3}.spm.stats.con.consess{30}.tcon.weights = [0 0 0 0 0 0 0 0 0 0 1 0 1 0 0 0 0 0];  
        matlabbatch{3}.spm.stats.con.consess{30}.tcon.sessrep = 'none';
        matlabbatch{3}.spm.stats.con.delete = 0;
        
        spm('defaults', 'FMRI');
        spm_jobman('run', matlabbatch);
        fileIn = fullfile(fmrifolderin,'Res*.nii'); system(['rm ' fileIn]);
        
        % Compute FD, DVARS, and motion representative values to reject outlier
        % subjects
        R(:,4:6) = R(:,4:6)*180/pi;
        Motion(id_s).MaxDisp = max(abs(R),[],1);
        Motion(id_s).MaxRelDisp = max(abs(diff(R,1,1)),[],1);
        Motion(id_s).MeanRMS = rms(R,1);
        Motion(id_s).id = subjname;
        clear R path_vo* Ons* volumes*
        save([fmrifolderin '/All_motionparam_StopTask.mat'],'Motion');
        clear matlabbatch;
        
    end
end