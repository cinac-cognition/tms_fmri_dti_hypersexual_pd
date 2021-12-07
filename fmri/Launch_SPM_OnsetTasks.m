%fMRI behavioural onsets and volumes integration
clear; close all; clc;

addpath(genpath([path 'spm12']));

source_dir = 'raw_data_path';
struct_dir = 'T1_data_path';
res_dir =    'spm_data_results';

for id_dataset = 1:5
    switch id_dataset
        case 1
            source_dir = '/data_path/ICD/OFF';
            struct_dir = '/data_path/ICD/Structural';
            res_dir =    '/data_path/ICD/ICDOFF';
            subjs = dir([source_dir '/01*']);
        case 2
            source_dir = '/data_path/ICD/ON';
            struct_dir = '/data_path/ICD/Structural';
            res_dir =    '/data_path/ICD/ICDON';
            subjs = dir([source_dir '/01*']);
        case 3
            source_dir = '/data_path/PD/OFF';
            struct_dir = '/data_path/PD/Structural';
            res_dir =    '/data_path/PD/PDOFF';
            subjs = dir([source_dir '/02*']);
        case 4
            source_dir = '/data_path/PD/ON';
            struct_dir = '/data_path/PD/Structural';
            res_dir =    '/data_path/PD/PDON';
            subjs = dir([source_dir '/02*']);
        case 5
            source_dir = '/data_path/controls';
            struct_dir = '/data_path/controls/Structural';
            res_dir =    '/data_path/controls1';
            subjs = dir([source_dir '/03*']);
    end
    
    if ~exist(res_dir,'dir'); mkdir(res_dir); end
    
    if 1
        for id_sub = 1:length(subjs)
            
            subjsname = subjs(id_sub).name;
            fprintf('Doing subject %s......\n',subjsname);
            subject_dir = fullfile(source_dir,subjsname);
            
            EPIfolders = dir([subject_dir '/EPI_BLOQUE_0*']);
            if length(EPIfolders)>5
                error('too many blocks for this subject...');
            end
            res_dir_subj = fullfile(res_dir,subjsname);
            if ~exist(res_dir_subj,'dir'); mkdir(res_dir_subj); end
            
            %--------------------------------------------------------------------------
            % Prepare files and folders
            fMRIdir_all = [res_dir_subj '/fMRI_all'];
            mkdir(fMRIdir_all);
            cont = 0;
            for i=1:length(EPIfolders)
                sub_tasks{i} = ['/fMRI_B' num2str(i)]; % Store each block folder
                fMRIdir = fullfile(res_dir_subj,sub_tasks{i});
                if exist(fMRIdir,'dir')
                    system(['rm -Rf ' fMRIdir]);
                    mkdir(fMRIdir);
                else
                    mkdir(fMRIdir);
                end
                
                % Convert dicom to nifti
                runline = sprintf('dcm2nii -o %s %s',fMRIdir,[subject_dir '/' EPIfolders(i).name]);
                [status,res] = system(runline);
                
                % Split 4D volumes into 3Ds
                filein = dir([fMRIdir '/EPI*']);
                runline = sprintf('fslsplit %s %s -t',fullfile(fMRIdir,filein.name),[fMRIdir '/vol']);
                [status,res] = system(runline); delete(fullfile(fMRIdir,filein.name));
                system(['gunzip ' fMRIdir '/*nii.gz']);
                vols = dir([fMRIdir '/vol*']);
                
                % Store number of volumes per block
                Nvols_x_bloq{id_sub}(i,1) = length(vols);
                % Store acquisition time per block related to the TR (2000ms)
                Tvols_x_bloq{id_sub}(i,1) = 2000*length(vols);
                
                for id_v = 1:length(vols)
                    cont = cont + 1;
                    system(['mv ' fMRIdir '/vol' num2str(id_v-1,'%04d') '.nii ' fMRIdir_all '/vol_T' num2str(i) '_' num2str(cont,'%04d') '.nii']);
                end
                system(['rm -Rf ' fMRIdir]);
            end
            %------------------------------------------------------------------
            % Read triggers file
            triggers_file = dir([subject_dir '/StopSignal*.txt']);
            triggers = importdata([subject_dir '/' triggers_file.name],'\t',1);
            %------------------------------------------------------------------
            % Fix Onset of the triggers file where task begins
            cont_B = 0;
            blockvector = [];
            for id_t = 1:size(triggers.data,1)
                sample = triggers.textdata{id_t+1,1};
                if strcmp(sample,'1')
                    cont_B = cont_B + 1;
                    blockvector(id_t) = cont_B;
                    trigger_init(cont_B) = triggers.data(id_t,6); % Fixation Onset - Stores per block the starting point
                    trigger_rowinit(cont_B) = id_t; % Stores per block the first row within the triggers file
                    if cont_B > 1
                        trigger_finish(cont_B-1) = triggers.data(id_t-1,12) + 1000; % ITIOnset + DurITI
                        trigger_rowfinish(cont_B-1) = id_t-1; % Stores per block the last row within the triggers file
                    end
                else
                    blockvector(id_t) = cont_B;
                end
            end
            trigger_finish(cont_B) = triggers.data(id_t,12) + 1000;
            trigger_rowfinish(cont_B) = id_t;
            trigger_flag_reset = [0 trigger_init(2:end)] < [0 trigger_finish(1:end-1)];
            time_x_block = trigger_finish - trigger_init;
            excessTime_x_block = Tvols_x_bloq{id_sub}' - time_x_block;
            excessVols_x_block = floor(abs(excessTime_x_block./2000)).*sign(excessTime_x_block);
            
            % Delete excess volumes
            for i_Bloque = 1:length(EPIfolders)
                if excessVols_x_block(i_Bloque) > 0
                    vols = dir([fMRIdir_all '/vol_T' num2str(i_Bloque) '*']);
                    for i_del = 1:excessVols_x_block(i_Bloque)
                        delete([fMRIdir_all '/' vols(end-i_del+1).name]);
                    end
                end
            end
            
            alltriggers2remove = [];
            acquisitionoffset = [];
            for i_Bloque = 1:length(EPIfolders)
                if excessTime_x_block(i_Bloque) < 0
                    trigger_finish(i_Bloque) = trigger_init(i_Bloque) + Tvols_x_bloq{id_sub}(i_Bloque);
                    Tfinish = nan(1,(1+trigger_rowfinish(i_Bloque)-trigger_rowinit(i_Bloque)));
                    for id_row = 1:(1+trigger_rowfinish(i_Bloque)-trigger_rowinit(i_Bloque))
                        Tfinish(id_row) = triggers.data(trigger_rowinit(i_Bloque)+id_row-1,12) + 1000;
                    end
                    triggers2remove(i_Bloque) = find(trigger_finish(i_Bloque)<Tfinish,1,'first');
                    trigger_finish_trial(i_Bloque) = Tfinish(triggers2remove(i_Bloque)-1);
                    alltriggers2remove = [alltriggers2remove (trigger_rowinit(i_Bloque)+triggers2remove(i_Bloque)-1):trigger_rowfinish(i_Bloque)];
                    acquisitionoffset(i_Bloque) = trigger_init(i_Bloque) + Tvols_x_bloq{id_sub}(i_Bloque) - trigger_finish_trial(i_Bloque);
                    
                else
                    triggers2remove(i_Bloque) = 0;
                end
            end
            clear Tfinish;
            %------------------------------------------------------------------
            triggers.data(alltriggers2remove,:) = [];
            triggers.textdata(alltriggers2remove,:) = [];
            blockvector(alltriggers2remove) = [];
            delta_finishinit = [0 trigger_init(2:end) - trigger_finish(1:end-1)];
            interblock_times = trigger_init(1) + [cumsum(delta_finishinit)];
            offset_reset = zeros(1,length(trigger_flag_reset));
            for id_reset = 1:length(EPIfolders)
                if trigger_flag_reset(id_reset) == 1
                    offset_reset(id_reset:end) = offset_reset(id_reset:end) ...
                        - (trigger_init(id_reset)) + trigger_finish(id_reset-1);
                end
            end
            delta_finishinit(trigger_flag_reset) = 1;
            interblock_times = trigger_init(1) + [cumsum(delta_finishinit)] - offset_reset;
            interblock_times_vec = [];
            Trialnames = {'Go' 'Go2' 'SSD1' 'SSD2' 'SSD3' 'SSD4' 'SSDa' 'SSDb' 'SSDc' 'SSDd'};
            PhotoCOND = triggers.data(:,1);
            cont_go = 0; cont_go2 = 0;
            cont_ssd1 = 0; cont_ssd2 = 0; cont_ssd3 = 0; cont_ssd4 = 0;
            cont_ssda = 0; cont_ssdb = 0; cont_ssdc = 0; cont_ssdd = 0;
            cont_ssd1_f = 0; cont_ssd2_f = 0; cont_ssd3_f = 0; cont_ssd4_f = 0;
            cont_ssda_f = 0; cont_ssdb_f = 0; cont_ssdc_f = 0; cont_ssdd_f = 0;
            pic_ero = 0; pic_noero = 0;
            StopOnset_fail = cell(10,1);
            StopOnset = cell(10,1);
            %-------------------------------------------
            for id_t = 1:size(triggers.data,1)
                offset = interblock_times(blockvector(id_t));
                condition = triggers.textdata{id_t+1,2};
                if strcmp(condition,'Go2')
                    if PhotoCOND(id_t) == 1
                        condition = 'Go';
                    end
                end
                if strcmp(condition,'Go')
                    if PhotoCOND(id_t) == 2
                        condition = 'Go2';
                    end
                end
                FixOnset(id_t,1) = (triggers.data(id_t,6) - offset)/1000;
                FixDuration(id_t,1) = (triggers.data(id_t,8) + 3000)/1000;
                % If it was accurate
                if (triggers.data(id_t,4)) == 1
                    switch condition
                        case 'Go' %Go erotic
                            pic_ero = pic_ero + 1;
                            cont_go = cont_go + 1;
                            trialcond = 1;
                            PicOnset_ero(pic_ero,1) = (triggers.data(id_t,9) - offset)/1000;
                            GoOnset_ero(cont_go,1) = (triggers.data(id_t,10) - offset)/1000;
                        case 'Go2' %Go non-erotic
                            cont_go2 = cont_go2 + 1;
                            pic_noero = pic_noero + 1;
                            trialcond = 2;
                            PicOnset_noero(pic_noero,1) = (triggers.data(id_t,9) - offset)/1000;
                            GoOnset_noero(cont_go2,1) = (triggers.data(id_t,10) - offset)/1000;
                        case 'SSD1' %SSD erotic
                            cont_ssd1 = cont_ssd1 + 1;
                            pic_ero = pic_ero + 1;
                            trialcond = 3;
                            PicOnset_ero(pic_ero,1) = (triggers.data(id_t,9) - offset)/1000;
                            StopOnset{trialcond,1}(cont_ssd1,1) = (triggers.data(id_t,11) - offset)/1000;
                        case 'SSD2' %SSD erotic
                            cont_ssd2 = cont_ssd2 + 1;
                            pic_ero = pic_ero + 1;
                            trialcond = 4;
                            PicOnset_ero(pic_ero,1) = (triggers.data(id_t,9) - offset)/1000;
                            StopOnset{trialcond,1}(cont_ssd2,1) = (triggers.data(id_t,11) - offset)/1000;
                        case 'SSD3' %SSD erotic
                            cont_ssd3 = cont_ssd3 + 1;
                            pic_ero = pic_ero + 1;
                            trialcond = 5;
                            PicOnset_ero(pic_ero,1) = (triggers.data(id_t,9) - offset)/1000;
                            StopOnset{trialcond,1}(cont_ssd3,1) = (triggers.data(id_t,11) - offset)/1000;
                        case 'SSD4' %SSD erotic
                            cont_ssd4 = cont_ssd4 + 1;
                            pic_ero = pic_ero + 1;
                            trialcond = 6;
                            PicOnset_ero(pic_ero,1) = (triggers.data(id_t,9) - offset)/1000;
                            StopOnset{trialcond,1}(cont_ssd4,1) = (triggers.data(id_t,11) - offset)/1000;
                        case 'SSDa' %SSD non-erotic
                            cont_ssda = cont_ssda + 1;
                            pic_noero = pic_noero + 1;
                            trialcond = 7;
                            PicOnset_noero(pic_noero,1) = (triggers.data(id_t,9) - offset)/1000;
                            StopOnset{trialcond,1}(cont_ssda,1) = (triggers.data(id_t,11) - offset)/1000;
                        case 'SSDb' %SSD non-erotic
                            cont_ssdb = cont_ssdb + 1;
                            pic_noero = pic_noero + 1;
                            trialcond = 8;
                            PicOnset_noero(pic_noero,1) = (triggers.data(id_t,9) - offset)/1000;
                            StopOnset{trialcond,1}(cont_ssdb,1) = (triggers.data(id_t,11) - offset)/1000;
                        case 'SSDc' %SSD non-erotic
                            cont_ssdc = cont_ssdc + 1;
                            pic_noero = pic_noero + 1;
                            trialcond = 9;
                            PicOnset_noero(pic_noero,1) = (triggers.data(id_t,9) - offset)/1000;
                            StopOnset{trialcond,1}(cont_ssdc,1) = (triggers.data(id_t,11) - offset)/1000;
                        case 'SSDd' %SSD non-erotic
                            cont_ssdd = cont_ssdd + 1;
                            pic_noero = pic_noero + 1;
                            trialcond = 10;
                            PicOnset_noero(pic_noero,1) = (triggers.data(id_t,9) - offset)/1000;
                            StopOnset{trialcond,1}(cont_ssdd,1) = (triggers.data(id_t,11) - offset)/1000;
                        otherwise
                            error('Trial Conditions Not Recognized');
                    end
                end
            end
            
            StopOnsetOk_ero = sort([StopOnset{3};StopOnset{4};StopOnset{5};StopOnset{6}],'ascend');
            StopOnsetOk_noero = sort([StopOnset{7};StopOnset{8};StopOnset{9};StopOnset{10}],'ascend');
            StopOnsetFail_ero = sort([StopOnset_fail{3};StopOnset_fail{4};StopOnset_fail{5};StopOnset_fail{6}],'ascend');
            StopOnsetFail_noero = sort([StopOnset_fail{7};StopOnset_fail{8};StopOnset_fail{9};StopOnset_fail{10}],'ascend');
            
            StopOnsetOk_ero(isnan(StopOnsetOk_ero)) = [];
            StopOnsetOk_noero(isnan(StopOnsetOk_noero)) = [];
            StopOnsetFail_ero(isnan(StopOnsetFail_ero)) = [];
            StopOnsetFail_noero(isnan(StopOnsetFail_noero)) = [];
            
            PicOnset_ero(isnan(PicOnset_ero)) = [];
            PicOnset_noero(isnan(PicOnset_noero)) = [];            
            GoOnset_ero(isnan(GoOnset_ero)) = [];
            GoOnset_noero(isnan(GoOnset_noero)) = [];            
            FixOnset(isnan(FixOnset)) = [];
            FixDuration(isnan(FixDuration)) = [];            
            save(fullfile(fMRIdir_all,'Onsets.mat'),'StopOnsetOk_ero','StopOnsetOk_noero','StopOnsetFail_ero','StopOnsetFail_noero','PicOnset_ero','PicOnset_noero','GoOnset_ero','GoOnset_noero','StopOnset','StopOnset_fail','FixOnset','FixDuration');
            save(fullfile(fMRIdir_all,'SSD_Onsets.mat'),'StopOnsetOk_ero','StopOnsetOk_noero');
        end
    end    
end

