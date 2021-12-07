%% Creation of DCM models and BMS

clear all; close all; clc;

root{1} = 'subject_data';
cd(root)
spm('defaults','FMRI');spm_jobman('initcfg');
ana_dir='fMRI_all_Stat';
clear matlabbatch;
matlabbatch=[];

for i_fold = 1
    root_pat = root{i_fold};
    sb = dir([root_pat '/0*']);
    
    for j=1:length(sb)
        output_dir = fullfile(root{i_fold},sb(j).name,ana_dir,'DCM_models');
        if ~exist(output_dir,'dir'); mkdir(output_dir); end
        %% Full model; comments explain each DCM matrix only for this model as applies to the others
        clear DCM
        load(fullfile(root_pat,sb(j).name,ana_dir,'SPM.mat'));
        % Regions of interest
        load(fullfile(root_pat,sb(j).name,ana_dir,'VOI_VTA.mat'),'xY');
        DCM.xY(1) = xY;
        load(fullfile(root_pat,sb(j).name,ana_dir,'VOI_caudate.mat'),'xY');
        DCM.xY(2) = xY;
        load(fullfile(root_pat,sb(j).name,ana_dir,'VOI_ACC.mat'),'xY');
        DCM.xY(3) = xY;
        load(fullfile(root_pat,sb(j).name,ana_dir,'VOI_preSMA.mat'),'xY');
        DCM.xY(4) = xY;
        DCM.n = length(DCM.xY);      % number of regions
        DCM.v = length(DCM.xY(1).u); % number of time points
        % Time series
        DCM.Y.dt  = SPM.xY.RT;
        DCM.Y.X0  = DCM.xY(1).X0;
        for i = 1:DCM.n
            DCM.Y.y(:,i)  = DCM.xY(i).u;
            DCM.Y.name{i} = DCM.xY(i).name;
        end
        DCM.Y.Q = spm_Ce(ones(1,DCM.n)*DCM.v);
        % Experimental inputs
        DCM.U.dt   =  SPM.Sess.U(2).dt;
        DCM.U.name = cellstr([SPM.Sess.U(6).name, SPM.Sess.U(2).name, SPM.Sess.U(3).name]);
        DCM.U.u    = [SPM.Sess.U(6).u(1:end,1),SPM.Sess.U(2).u(1:end,1),SPM.Sess.U(3).u(1:end,1)];
        
        % DCM parameters and options
        DCM.delays = repmat(SPM.xY.RT/2,DCM.n,1);
        DCM.TE     = 0.04;
        DCM.options.nonlinear  = 0;
        DCM.options.two_state  = 1;
        DCM.options.stochastic = 0;
        DCM.options.nograph    = 1;
        % Connectivity matrices: VTA,Caudate,ACC,preSMA
        % Create fixed connections DCM.a
        DCM.a=[1 1 1 1;
            1 1 1 1;
            1 1 1 1;
            1 1 1 1];
        % B matrix: task variables on connections
        DCM.b = zeros(3,3,3);
        %1.Erotic Image connections
        DCM.b(:,:,1) = [1 1 1 1;
            1 1 1 1;
            1 1 1 1;
            1 1 1 1];
        %2.stop task erotic connections
        DCM.b(:,:,2) = [1 1 1 1;
            1 1 1 1;
            1 1 1 1;
            1 1 1 1];
        %3.stop task Non-erotic connections
        DCM.b(:,:,3) = [1 1 1 1;
            1 1 1 1;
            1 1 1 1;
            1 1 1 1];
        % Task-effects over each region:
        % rows: VTA,Caudate,ACC,preSMA as rows
        % columns: 1st erotic image, 2nd stop erotic, 3rd stop non-ero
        DCM.c =[1 1 1 1;
            1 1 1 1;
            1 1 1 1;
            1 1 1 1];

        save(fullfile(output_dir,'DCM_full.mat'),'DCM');
        % DCM Estimation
        clear matlabbatch
        matlabbatch{1}.spm.dcm.fmri.estimate.dcmmat = {fullfile(output_dir,'DCM_full.mat')};
        spm_jobman('run',matlabbatch);
        
        %% model1
        clear DCM
        load(fullfile(root_pat,sb(j).name,ana_dir,'SPM.mat'));
        % Regions of interest
        load(fullfile(root_pat,sb(j).name,ana_dir,'VOI_VTA.mat'),'xY');
        DCM.xY(1) = xY;
        load(fullfile(root_pat,sb(j).name,ana_dir,'VOI_caudate.mat'),'xY');
        DCM.xY(2) = xY;
        load(fullfile(root_pat,sb(j).name,ana_dir,'VOI_ACC.mat'),'xY');
        DCM.xY(3) = xY;
        load(fullfile(root_pat,sb(j).name,ana_dir,'VOI_preSMA.mat'),'xY');
        DCM.xY(4) = xY;
        DCM.n = length(DCM.xY);      % number of regions
        DCM.v = length(DCM.xY(1).u); % number of time points
        % Time series
        DCM.Y.dt  = SPM.xY.RT;
        DCM.Y.X0  = DCM.xY(1).X0;
        for i = 1:DCM.n
            DCM.Y.y(:,i)  = DCM.xY(i).u;
            DCM.Y.name{i} = DCM.xY(i).name;
        end
        DCM.Y.Q = spm_Ce(ones(1,DCM.n)*DCM.v);
        % Experimental inputs
        DCM.U.dt   =  SPM.Sess.U(2).dt; 
        DCM.U.name = cellstr([SPM.Sess.U(6).name, SPM.Sess.U(2).name, SPM.Sess.U(3).name]);
        DCM.U.u    = [SPM.Sess.U(6).u(1:end,1),SPM.Sess.U(2).u(1:end,1),SPM.Sess.U(3).u(1:end,1)];
        % DCM parameters and options
        DCM.delays = repmat(SPM.xY.RT/2,DCM.n,1);
        DCM.TE     = 0.04;
        DCM.options.nonlinear  = 0;
        DCM.options.two_state  = 1;
        DCM.options.stochastic = 0;
        DCM.options.nograph    = 1;
        % Connectivity matrices
        DCM.a=[1 0 0 0;
            1 1 0 0;
            1 0 1 0;
            0 1 1 1];
        DCM.b = zeros(4,4,3);
        DCM.b(4,3,1) = 1;DCM.b(:,:,2) = 0;DCM.b(:,:,3) = 0;
        DCM.c =[0 0 0;
            0 0 0;
            1 0 0;
            0 1 1];
        save(fullfile(output_dir,'DCM_1.mat'),'DCM');
        % % DCM Estimation
        clear matlabbatch
        matlabbatch{1}.spm.dcm.fmri.estimate.dcmmat = {fullfile(output_dir,'DCM_1.mat')};
        spm_jobman('run',matlabbatch);
        %% model2
        clear DCM
        load(fullfile(root_pat,sb(j).name,ana_dir,'SPM.mat'));
        % Regions of interest
        load(fullfile(root_pat,sb(j).name,ana_dir,'VOI_VTA.mat'),'xY');
        DCM.xY(1) = xY;
        load(fullfile(root_pat,sb(j).name,ana_dir,'VOI_caudate.mat'),'xY');
        DCM.xY(2) = xY;
        load(fullfile(root_pat,sb(j).name,ana_dir,'VOI_ACC.mat'),'xY');
        DCM.xY(3) = xY;
        load(fullfile(root_pat,sb(j).name,ana_dir,'VOI_preSMA.mat'),'xY');
        DCM.xY(4) = xY;
        DCM.n = length(DCM.xY);      % number of regions
        DCM.v = length(DCM.xY(1).u); % number of time points
        % Time series
        DCM.Y.dt  = SPM.xY.RT;
        DCM.Y.X0  = DCM.xY(1).X0;
        for i = 1:DCM.n
            DCM.Y.y(:,i)  = DCM.xY(i).u;
            DCM.Y.name{i} = DCM.xY(i).name;
        end
        DCM.Y.Q = spm_Ce(ones(1,DCM.n)*DCM.v);
        % Experimental inputs
        DCM.U.dt   =  SPM.Sess.U(2).dt; 
        DCM.U.name = cellstr([SPM.Sess.U(6).name, SPM.Sess.U(2).name, SPM.Sess.U(3).name]);
        DCM.U.u    = [SPM.Sess.U(6).u(1:end,1),SPM.Sess.U(2).u(1:end,1),SPM.Sess.U(3).u(1:end,1)];
        % DCM parameters and options
        DCM.delays = repmat(SPM.xY.RT/2,DCM.n,1);
        DCM.TE     = 0.04;
        DCM.options.nonlinear  = 0;
        DCM.options.two_state  = 1;
        DCM.options.stochastic = 0;
        DCM.options.nograph    = 1;
        % Connectivity matrices
        DCM.a=[1 0 0 0;
            0 1 0 0;
            0 0 1 0;
            0 1 1 1];
        DCM.b = zeros(4,4,3);
        DCM.b(4,3,1) = 1;DCM.b(:,:,2) = 0;DCM.b(:,:,3) = 0;
        DCM.c =[0 0 0;
            0 0 0;
            1 0 0;
            0 1 1];
        save(fullfile(output_dir,'DCM_2.mat'),'DCM');
        % % DCM Estimation
        clear matlabbatch
        matlabbatch{1}.spm.dcm.fmri.estimate.dcmmat = {fullfile(output_dir,'DCM_2.mat')};
        spm_jobman('run',matlabbatch);
        %% model3
        clear DCM
        load(fullfile(root_pat,sb(j).name,ana_dir,'SPM.mat'));
        % Regions of interest
        load(fullfile(root_pat,sb(j).name,ana_dir,'VOI_VTA.mat'),'xY');
        DCM.xY(1) = xY;
        load(fullfile(root_pat,sb(j).name,ana_dir,'VOI_caudate.mat'),'xY');
        DCM.xY(2) = xY;
        load(fullfile(root_pat,sb(j).name,ana_dir,'VOI_ACC.mat'),'xY');
        DCM.xY(3) = xY;
        load(fullfile(root_pat,sb(j).name,ana_dir,'VOI_preSMA.mat'),'xY');
        DCM.xY(4) = xY;
        DCM.n = length(DCM.xY);      % number of regions
        DCM.v = length(DCM.xY(1).u); % number of time points
        % Time series
        DCM.Y.dt  = SPM.xY.RT;
        DCM.Y.X0  = DCM.xY(1).X0;
        for i = 1:DCM.n
            DCM.Y.y(:,i)  = DCM.xY(i).u;
            DCM.Y.name{i} = DCM.xY(i).name;
        end
        DCM.Y.Q = spm_Ce(ones(1,DCM.n)*DCM.v);
        % Experimental inputs
        DCM.U.dt   =  SPM.Sess.U(2).dt;
        DCM.U.name = cellstr([SPM.Sess.U(6).name, SPM.Sess.U(2).name, SPM.Sess.U(3).name]);
        DCM.U.u    = [SPM.Sess.U(6).u(1:end,1),SPM.Sess.U(2).u(1:end,1),SPM.Sess.U(3).u(1:end,1)];
        % DCM parameters and options
        DCM.delays = repmat(SPM.xY.RT/2,DCM.n,1);
        DCM.TE     = 0.04;
        DCM.options.nonlinear  = 0;
        DCM.options.two_state  = 1;
        DCM.options.stochastic = 0;
        DCM.options.nograph    = 1;
        % Connectivity matrices
        DCM.a=[1 0 0 0;
            0 1 0 1;
            0 0 1 0;
            0 1 1 1];
        DCM.b = zeros(4,4,3);
        DCM.b(4,3,1) = 1;DCM.b(:,:,2) = 0;DCM.b(:,:,3) = 0;
        DCM.c =[0 0 0;
            0 0 0;
            1 0 0;
            0 1 1];
        save(fullfile(output_dir,'DCM_3.mat'),'DCM');
        % % DCM Estimation
        clear matlabbatch
        matlabbatch{1}.spm.dcm.fmri.estimate.dcmmat = {fullfile(output_dir,'DCM_3.mat')};
        spm_jobman('run',matlabbatch);
        
        %% model4
        clear DCM
        load(fullfile(root_pat,sb(j).name,ana_dir,'SPM.mat'));
        % Regions of interest
        load(fullfile(root_pat,sb(j).name,ana_dir,'VOI_VTA.mat'),'xY');
        DCM.xY(1) = xY;
        load(fullfile(root_pat,sb(j).name,ana_dir,'VOI_caudate.mat'),'xY');
        DCM.xY(2) = xY;
        load(fullfile(root_pat,sb(j).name,ana_dir,'VOI_ACC.mat'),'xY');
        DCM.xY(3) = xY;
        load(fullfile(root_pat,sb(j).name,ana_dir,'VOI_preSMA.mat'),'xY');
        DCM.xY(4) = xY;
        DCM.n = length(DCM.xY);      % number of regions
        DCM.v = length(DCM.xY(1).u); % number of time points
        % Time series
        DCM.Y.dt  = SPM.xY.RT;
        DCM.Y.X0  = DCM.xY(1).X0;
        for i = 1:DCM.n
            DCM.Y.y(:,i)  = DCM.xY(i).u;
            DCM.Y.name{i} = DCM.xY(i).name;
        end
        DCM.Y.Q = spm_Ce(ones(1,DCM.n)*DCM.v);
        % Experimental inputs
        DCM.U.dt   =  SPM.Sess.U(2).dt; 
        DCM.U.name = cellstr([SPM.Sess.U(6).name, SPM.Sess.U(2).name, SPM.Sess.U(3).name]);
        DCM.U.u    = [SPM.Sess.U(6).u(1:end,1),SPM.Sess.U(2).u(1:end,1),SPM.Sess.U(3).u(1:end,1)];
        % DCM parameters and options
        DCM.delays = repmat(SPM.xY.RT/2,DCM.n,1);
        DCM.TE     = 0.04;
        DCM.options.nonlinear  = 0;
        DCM.options.two_state  = 1;
        DCM.options.stochastic = 0;
        DCM.options.nograph    = 1;
        % Connectivity matrices
        DCM.a=[1 0 0 0;
            0 1 0 1;
            0 0 1 0;
            0 1 1 1];
        DCM.b = zeros(4,4,3);
        DCM.b(:,:,1) = 0;DCM.b(2,4,2) = 1;DCM.b(2,4,3) = 1;
        DCM.c =[0 0 0;
            0 0 0;
            1 0 0;
            0 1 1];
        save(fullfile(output_dir,'DCM_4.mat'),'DCM');
        % % DCM Estimation
        clear matlabbatch
        matlabbatch{1}.spm.dcm.fmri.estimate.dcmmat = {fullfile(output_dir,'DCM_4.mat')};
        spm_jobman('run',matlabbatch);
        
        %% model5
        clear DCM
        load(fullfile(root_pat,sb(j).name,ana_dir,'SPM.mat'));
        % Regions of interest
        load(fullfile(root_pat,sb(j).name,ana_dir,'VOI_VTA.mat'),'xY');
        DCM.xY(1) = xY;
        load(fullfile(root_pat,sb(j).name,ana_dir,'VOI_caudate.mat'),'xY');
        DCM.xY(2) = xY;
        load(fullfile(root_pat,sb(j).name,ana_dir,'VOI_ACC.mat'),'xY');
        DCM.xY(3) = xY;
        load(fullfile(root_pat,sb(j).name,ana_dir,'VOI_preSMA.mat'),'xY');
        DCM.xY(4) = xY;
        DCM.n = length(DCM.xY);      % number of regions
        DCM.v = length(DCM.xY(1).u); % number of time points
        % Time series
        DCM.Y.dt  = SPM.xY.RT;
        DCM.Y.X0  = DCM.xY(1).X0;
        for i = 1:DCM.n
            DCM.Y.y(:,i)  = DCM.xY(i).u;
            DCM.Y.name{i} = DCM.xY(i).name;
        end
        DCM.Y.Q = spm_Ce(ones(1,DCM.n)*DCM.v);
        % Experimental inputs
        DCM.U.dt   =  SPM.Sess.U(2).dt;
        DCM.U.name = cellstr([SPM.Sess.U(6).name, SPM.Sess.U(2).name, SPM.Sess.U(3).name]);
        DCM.U.u    = [SPM.Sess.U(6).u(1:end,1),SPM.Sess.U(2).u(1:end,1),SPM.Sess.U(3).u(1:end,1)];
        % DCM parameters and options
        DCM.delays = repmat(SPM.xY.RT/2,DCM.n,1);
        DCM.TE     = 0.04;
        DCM.options.nonlinear  = 0;
        DCM.options.two_state  = 1;
        DCM.options.stochastic = 0;
        DCM.options.nograph    = 1;
        % Connectivity matrices
        DCM.a=[1 0 0 0;
            0 1 0 1;
            0 1 1 0;
            0 1 1 1];
        DCM.b = zeros(4,4,3);
        DCM.b(4,3,1) = 1;DCM.b(:,:,2) = 0;DCM.b(:,:,3) = 0;
        DCM.c =[0 0 0;
            0 0 0;
            1 0 0;
            0 1 1];
        % % DCM Estimation
        clear matlabbatch
        matlabbatch{1}.spm.dcm.fmri.estimate.dcmmat = {fullfile(output_dir,'DCM_5.mat')};
        spm_jobman('run',matlabbatch);
        %% model6
        clear DCM
        load(fullfile(root_pat,sb(j).name,ana_dir,'SPM.mat'));
        % Regions of interest
        load(fullfile(root_pat,sb(j).name,ana_dir,'VOI_VTA.mat'),'xY');
        DCM.xY(1) = xY;
        load(fullfile(root_pat,sb(j).name,ana_dir,'VOI_caudate.mat'),'xY');
        DCM.xY(2) = xY;
        load(fullfile(root_pat,sb(j).name,ana_dir,'VOI_ACC.mat'),'xY');
        DCM.xY(3) = xY;
        load(fullfile(root_pat,sb(j).name,ana_dir,'VOI_preSMA.mat'),'xY');
        DCM.xY(4) = xY;
        DCM.n = length(DCM.xY);      % number of regions
        DCM.v = length(DCM.xY(1).u); % number of time points
        % Time series
        DCM.Y.dt  = SPM.xY.RT;
        DCM.Y.X0  = DCM.xY(1).X0;
        for i = 1:DCM.n
            DCM.Y.y(:,i)  = DCM.xY(i).u;
            DCM.Y.name{i} = DCM.xY(i).name;
        end
        DCM.Y.Q = spm_Ce(ones(1,DCM.n)*DCM.v);
        % Experimental inputs
        DCM.U.dt   =  SPM.Sess.U(2).dt;
        DCM.U.name = cellstr([SPM.Sess.U(6).name, SPM.Sess.U(2).name, SPM.Sess.U(3).name]);
        DCM.U.u    = [SPM.Sess.U(6).u(1:end,1),SPM.Sess.U(2).u(1:end,1),SPM.Sess.U(3).u(1:end,1)];
        % DCM parameters and options
        DCM.delays = repmat(SPM.xY.RT/2,DCM.n,1);
        DCM.TE     = 0.04;
        DCM.options.nonlinear  = 0;
        DCM.options.two_state  = 1;
        DCM.options.stochastic = 0;
        DCM.options.nograph    = 1;
        % Connectivity matrices
        DCM.a=[1 0 0 0;
            0 1 0 1;
            0 1 1 0;
            0 1 1 1];
        DCM.b = zeros(4,4,3);
        DCM.b(:,:,1) = 0;DCM.b(2,4,2) =1;DCM.b(2,4,3) = 1;
        DCM.c =[0 0 0;
            0 0 0;
            1 0 0;
            0 1 1];
        save(fullfile(output_dir,'DCM_6.mat'),'DCM');
        % % DCM Estimation
        clear matlabbatch
        matlabbatch{1}.spm.dcm.fmri.estimate.dcmmat = {fullfile(output_dir,'DCM_6.mat')};
        spm_jobman('run',matlabbatch);
        %% model7
        clear DCM
        load(fullfile(root_pat,sb(j).name,ana_dir,'SPM.mat'));
        % Regions of interest
        load(fullfile(root_pat,sb(j).name,ana_dir,'VOI_VTA.mat'),'xY');
        DCM.xY(1) = xY;
        load(fullfile(root_pat,sb(j).name,ana_dir,'VOI_caudate.mat'),'xY');
        DCM.xY(2) = xY;
        load(fullfile(root_pat,sb(j).name,ana_dir,'VOI_ACC.mat'),'xY');
        DCM.xY(3) = xY;
        load(fullfile(root_pat,sb(j).name,ana_dir,'VOI_preSMA.mat'),'xY');
        DCM.xY(4) = xY;
        DCM.n = length(DCM.xY);      % number of regions
        DCM.v = length(DCM.xY(1).u); % number of time points
        % Time series
        DCM.Y.dt  = SPM.xY.RT;
        DCM.Y.X0  = DCM.xY(1).X0;
        for i = 1:DCM.n
            DCM.Y.y(:,i)  = DCM.xY(i).u;
            DCM.Y.name{i} = DCM.xY(i).name;
        end
        DCM.Y.Q = spm_Ce(ones(1,DCM.n)*DCM.v);
        % Experimental inputs
        DCM.U.dt   =  SPM.Sess.U(2).dt;
        DCM.U.name = cellstr([SPM.Sess.U(6).name, SPM.Sess.U(2).name, SPM.Sess.U(3).name]);
        DCM.U.u    = [SPM.Sess.U(6).u(1:end,1),SPM.Sess.U(2).u(1:end,1),SPM.Sess.U(3).u(1:end,1)];
        % DCM parameters and options
        DCM.delays = repmat(SPM.xY.RT/2,DCM.n,1);
        DCM.TE     = 0.04;
        DCM.options.nonlinear  = 0;
        DCM.options.two_state  = 1;
        DCM.options.stochastic = 0;
        DCM.options.nograph    = 1;
        % Connectivity matrices
        DCM.a=[1 0 0 0;
            1 1 0 0;
            1 1 1 0;
            0 1 1 1];
        DCM.b = zeros(4,4,3);
        DCM.b(4,3,1) = 1;DCM.b(:,:,2) = 0;DCM.b(:,:,3) = 0;
        DCM.c =[0 0 0;
            0 0 0;
            1 0 0;
            0 1 1];
        save(fullfile(output_dir,'DCM_7.mat'),'DCM');
        % % DCM Estimation
        clear matlabbatch
        matlabbatch{1}.spm.dcm.fmri.estimate.dcmmat = {fullfile(output_dir,'DCM_7.mat')};
        spm_jobman('run',matlabbatch);
        %% model8
        clear DCM
        load(fullfile(root_pat,sb(j).name,ana_dir,'SPM.mat'));
        % Regions of interest
        load(fullfile(root_pat,sb(j).name,ana_dir,'VOI_VTA.mat'),'xY');
        DCM.xY(1) = xY;
        load(fullfile(root_pat,sb(j).name,ana_dir,'VOI_caudate.mat'),'xY');
        DCM.xY(2) = xY;
        load(fullfile(root_pat,sb(j).name,ana_dir,'VOI_ACC.mat'),'xY');
        DCM.xY(3) = xY;
        load(fullfile(root_pat,sb(j).name,ana_dir,'VOI_preSMA.mat'),'xY');
        DCM.xY(4) = xY;
        DCM.n = length(DCM.xY);      % number of regions
        DCM.v = length(DCM.xY(1).u); % number of time points
        % Time series
        DCM.Y.dt  = SPM.xY.RT;
        DCM.Y.X0  = DCM.xY(1).X0;
        for i = 1:DCM.n
            DCM.Y.y(:,i)  = DCM.xY(i).u;
            DCM.Y.name{i} = DCM.xY(i).name;
        end
        DCM.Y.Q = spm_Ce(ones(1,DCM.n)*DCM.v);
        % Experimental inputs
        DCM.U.dt   =  SPM.Sess.U(2).dt;
        DCM.U.name = cellstr([SPM.Sess.U(6).name, SPM.Sess.U(2).name, SPM.Sess.U(3).name]);
        DCM.U.u    = [SPM.Sess.U(6).u(1:end,1),SPM.Sess.U(2).u(1:end,1),SPM.Sess.U(3).u(1:end,1)];
        % DCM parameters and options
        DCM.delays = repmat(SPM.xY.RT/2,DCM.n,1);
        DCM.TE     = 0.04;
        DCM.options.nonlinear  = 0;
        DCM.options.two_state  = 1;
        DCM.options.stochastic = 0;
        DCM.options.nograph    = 1;
        % Connectivity matrices
        DCM.a=[1 0 1 0;
            1 1 0 0;
            1 1 1 0;
            0 1 1 1];
        DCM.b = zeros(4,4,3);
        DCM.b(4,3,1) = 1;DCM.b(:,:,2) = 0;DCM.b(:,:,3) = 0;
        DCM.c =[0 0 0;
            0 0 0;
            1 0 0;
            0 1 1];
        save(fullfile(output_dir,'DCM_8.mat'),'DCM');
        % % DCM Estimation
        clear matlabbatch
        matlabbatch{1}.spm.dcm.fmri.estimate.dcmmat = {fullfile(output_dir,'DCM_8.mat')};
        spm_jobman('run',matlabbatch);
        %% model9
        clear DCM
        load(fullfile(root_pat,sb(j).name,ana_dir,'SPM.mat'));
        % Regions of interest
        load(fullfile(root_pat,sb(j).name,ana_dir,'VOI_VTA.mat'),'xY');
        DCM.xY(1) = xY;
        load(fullfile(root_pat,sb(j).name,ana_dir,'VOI_caudate.mat'),'xY');
        DCM.xY(2) = xY;
        load(fullfile(root_pat,sb(j).name,ana_dir,'VOI_ACC.mat'),'xY');
        DCM.xY(3) = xY;
        load(fullfile(root_pat,sb(j).name,ana_dir,'VOI_preSMA.mat'),'xY');
        DCM.xY(4) = xY;
        DCM.n = length(DCM.xY);      % number of regions
        DCM.v = length(DCM.xY(1).u); % number of time points
        % Time series
        DCM.Y.dt  = SPM.xY.RT;
        DCM.Y.X0  = DCM.xY(1).X0;
        for i = 1:DCM.n
            DCM.Y.y(:,i)  = DCM.xY(i).u;
            DCM.Y.name{i} = DCM.xY(i).name;
        end
        DCM.Y.Q = spm_Ce(ones(1,DCM.n)*DCM.v);
        % Experimental inputs
        DCM.U.dt   =  SPM.Sess.U(2).dt; 
        DCM.U.name = cellstr([SPM.Sess.U(6).name, SPM.Sess.U(2).name, SPM.Sess.U(3).name]);
        DCM.U.u    = [SPM.Sess.U(6).u(1:end,1),SPM.Sess.U(2).u(1:end,1),SPM.Sess.U(3).u(1:end,1)];
        % DCM parameters and options
        DCM.delays = repmat(SPM.xY.RT/2,DCM.n,1);
        DCM.TE     = 0.04;
        DCM.options.nonlinear  = 0;
        DCM.options.two_state  = 1;
        DCM.options.stochastic = 0;
        DCM.options.nograph    = 1;
        % Connectivity matrices
        DCM.a=[1 0 1 0;
            1 1 0 0;
            1 1 1 0;
            0 1 1 1];
        DCM.b = zeros(4,4,3);
        DCM.b(1,3,1) = 1;DCM.b(:,:,2) = 0;DCM.b(:,:,3) = 0;
        DCM.c =[0 0 0;
            0 0 0;
            1 0 0;
            0 1 1];
        save(fullfile(output_dir,'DCM_9.mat'),'DCM');
        % % DCM Estimation
        clear matlabbatch
        matlabbatch{1}.spm.dcm.fmri.estimate.dcmmat = {fullfile(output_dir,'DCM_9.mat')};
        spm_jobman('run',matlabbatch);
        %% model10
        clear DCM
        load(fullfile(root_pat,sb(j).name,ana_dir,'SPM.mat'));
        % Regions of interest
        load(fullfile(root_pat,sb(j).name,ana_dir,'VOI_VTA.mat'),'xY');
        DCM.xY(1) = xY;
        load(fullfile(root_pat,sb(j).name,ana_dir,'VOI_caudate.mat'),'xY');
        DCM.xY(2) = xY;
        load(fullfile(root_pat,sb(j).name,ana_dir,'VOI_ACC.mat'),'xY');
        DCM.xY(3) = xY;
        load(fullfile(root_pat,sb(j).name,ana_dir,'VOI_preSMA.mat'),'xY');
        DCM.xY(4) = xY;
        DCM.n = length(DCM.xY);      % number of regions
        DCM.v = length(DCM.xY(1).u); % number of time points
        % Time series
        DCM.Y.dt  = SPM.xY.RT;
        DCM.Y.X0  = DCM.xY(1).X0;
        for i = 1:DCM.n
            DCM.Y.y(:,i)  = DCM.xY(i).u;
            DCM.Y.name{i} = DCM.xY(i).name;
        end
        DCM.Y.Q = spm_Ce(ones(1,DCM.n)*DCM.v);
        % Experimental inputs
        DCM.U.dt   =  SPM.Sess.U(2).dt; 
        DCM.U.name = cellstr([SPM.Sess.U(6).name, SPM.Sess.U(2).name, SPM.Sess.U(3).name]);
        DCM.U.u    = [SPM.Sess.U(6).u(1:end,1),SPM.Sess.U(2).u(1:end,1),SPM.Sess.U(3).u(1:end,1)];
        % DCM parameters and options
        DCM.delays = repmat(SPM.xY.RT/2,DCM.n,1);
        DCM.TE     = 0.04;
        DCM.options.nonlinear  = 0;
        DCM.options.two_state  = 1;
        DCM.options.stochastic = 0;
        DCM.options.nograph    = 1;
        % Connectivity matrices
        DCM.a=[1 0 0 0;
            1 1 0 1;
            1 1 1 0;
            0 1 0 1];
        DCM.b = zeros(4,4,3);
        DCM.b(:,:,1) = 0;DCM.b(2,4,2) = 1;DCM.b(2,4,3) = 1;
        DCM.c =[0 0 0;
            0 0 0;
            1 0 0;
            0 1 1];
        save(fullfile(output_dir,'DCM_10.mat'),'DCM');
        % % DCM Estimation
        clear matlabbatch
        matlabbatch{1}.spm.dcm.fmri.estimate.dcmmat = {fullfile(output_dir,'DCM_10.mat')};
        spm_jobman('run',matlabbatch);
        %% model11
        clear DCM
        load(fullfile(root_pat,sb(j).name,ana_dir,'SPM.mat'));
        % Regions of interest
        load(fullfile(root_pat,sb(j).name,ana_dir,'VOI_VTA.mat'),'xY');
        DCM.xY(1) = xY;
        load(fullfile(root_pat,sb(j).name,ana_dir,'VOI_caudate.mat'),'xY');
        DCM.xY(2) = xY;
        load(fullfile(root_pat,sb(j).name,ana_dir,'VOI_ACC.mat'),'xY');
        DCM.xY(3) = xY;
        load(fullfile(root_pat,sb(j).name,ana_dir,'VOI_preSMA.mat'),'xY');
        DCM.xY(4) = xY;
        DCM.n = length(DCM.xY);      % number of regions
        DCM.v = length(DCM.xY(1).u); % number of time points
        % Time series
        DCM.Y.dt  = SPM.xY.RT;
        DCM.Y.X0  = DCM.xY(1).X0;
        for i = 1:DCM.n
            DCM.Y.y(:,i)  = DCM.xY(i).u;
            DCM.Y.name{i} = DCM.xY(i).name;
        end
        DCM.Y.Q = spm_Ce(ones(1,DCM.n)*DCM.v);
        % Experimental inputs
        DCM.U.dt   =  SPM.Sess.U(2).dt;
        DCM.U.name = cellstr([SPM.Sess.U(6).name, SPM.Sess.U(2).name, SPM.Sess.U(3).name]);
        DCM.U.u    = [SPM.Sess.U(6).u(1:end,1),SPM.Sess.U(2).u(1:end,1),SPM.Sess.U(3).u(1:end,1)];
        % DCM parameters and options
        DCM.delays = repmat(SPM.xY.RT/2,DCM.n,1);
        DCM.TE     = 0.04;
        DCM.options.nonlinear  = 0;
        DCM.options.two_state  = 1;
        DCM.options.stochastic = 0;
        DCM.options.nograph    = 1;
        % Connectivity matrices
        DCM.a=[1 0 0 0;
            1 1 0 0;
            1 0 1 1;
            0 1 1 1];
        DCM.b = zeros(4,4,3);
        DCM.b(4,3,1) = 1;DCM.b(:,:,2) = 0;DCM.b(:,:,3) = 0;
        DCM.c =[0 0 0;
            0 0 0;
            1 0 0;
            0 1 1];
        save(fullfile(output_dir,'DCM_11.mat'),'DCM');
        % % DCM Estimation
        clear matlabbatch
        matlabbatch{1}.spm.dcm.fmri.estimate.dcmmat = {fullfile(output_dir,'DCM_11.mat')};
        spm_jobman('run',matlabbatch);
        %% model12
        clear DCM
        load(fullfile(root_pat,sb(j).name,ana_dir,'SPM.mat'));
        % Regions of interest
        load(fullfile(root_pat,sb(j).name,ana_dir,'VOI_VTA.mat'),'xY');
        DCM.xY(1) = xY;
        load(fullfile(root_pat,sb(j).name,ana_dir,'VOI_caudate.mat'),'xY');
        DCM.xY(2) = xY;
        load(fullfile(root_pat,sb(j).name,ana_dir,'VOI_ACC.mat'),'xY');
        DCM.xY(3) = xY;
        load(fullfile(root_pat,sb(j).name,ana_dir,'VOI_preSMA.mat'),'xY');
        DCM.xY(4) = xY;
        DCM.n = length(DCM.xY);      % number of regions
        DCM.v = length(DCM.xY(1).u); % number of time points
        % Time series
        DCM.Y.dt  = SPM.xY.RT;
        DCM.Y.X0  = DCM.xY(1).X0;
        for i = 1:DCM.n
            DCM.Y.y(:,i)  = DCM.xY(i).u;
            DCM.Y.name{i} = DCM.xY(i).name;
        end
        DCM.Y.Q = spm_Ce(ones(1,DCM.n)*DCM.v);
        % Experimental inputs
        DCM.U.dt   =  SPM.Sess.U(2).dt;
        DCM.U.name = cellstr([SPM.Sess.U(6).name, SPM.Sess.U(2).name, SPM.Sess.U(3).name]);
        DCM.U.u    = [SPM.Sess.U(6).u(1:end,1),SPM.Sess.U(2).u(1:end,1),SPM.Sess.U(3).u(1:end,1)];
        
        % DCM parameters and options
        DCM.delays = repmat(SPM.xY.RT/2,DCM.n,1);
        DCM.TE     = 0.04;
        DCM.options.nonlinear  = 0;
        DCM.options.two_state  = 1;
        DCM.options.stochastic = 0;
        DCM.options.nograph    = 1;
        % Connectivity matrices
        DCM.a=[1 0 1 0;
            1 1 0 0;
            1 0 1 1;
            0 1 0 1];
        DCM.b = zeros(4,4,3);
        DCM.b(4,3,1) = 1;DCM.b(:,:,2) = 0;DCM.b(:,:,3) = 0;
        DCM.c =[0 0 0;
            0 0 0;
            1 0 0;
            0 1 1];
        save(fullfile(output_dir,'DCM_12.mat'),'DCM');
        % % DCM Estimation
        clear matlabbatch
        matlabbatch{1}.spm.dcm.fmri.estimate.dcmmat = {fullfile(output_dir,'DCM_12.mat')};
        spm_jobman('run',matlabbatch);
        %% model13
        clear DCM
        load(fullfile(root_pat,sb(j).name,ana_dir,'SPM.mat'));
        % Regions of interest
        load(fullfile(root_pat,sb(j).name,ana_dir,'VOI_VTA.mat'),'xY');
        DCM.xY(1) = xY;
        load(fullfile(root_pat,sb(j).name,ana_dir,'VOI_caudate.mat'),'xY');
        DCM.xY(2) = xY;
        load(fullfile(root_pat,sb(j).name,ana_dir,'VOI_ACC.mat'),'xY');
        DCM.xY(3) = xY;
        load(fullfile(root_pat,sb(j).name,ana_dir,'VOI_preSMA.mat'),'xY');
        DCM.xY(4) = xY;
        DCM.n = length(DCM.xY);      % number of regions
        DCM.v = length(DCM.xY(1).u); % number of time points
        % Time series
        DCM.Y.dt  = SPM.xY.RT;
        DCM.Y.X0  = DCM.xY(1).X0;
        for i = 1:DCM.n
            DCM.Y.y(:,i)  = DCM.xY(i).u;
            DCM.Y.name{i} = DCM.xY(i).name;
        end
        DCM.Y.Q = spm_Ce(ones(1,DCM.n)*DCM.v);
        % Experimental inputs
        DCM.U.dt   =  SPM.Sess.U(2).dt; 
        DCM.U.name = cellstr([SPM.Sess.U(6).name, SPM.Sess.U(2).name, SPM.Sess.U(3).name]);
        DCM.U.u    = [SPM.Sess.U(6).u(1:end,1),SPM.Sess.U(2).u(1:end,1),SPM.Sess.U(3).u(1:end,1)];
        % DCM parameters and options
        DCM.delays = repmat(SPM.xY.RT/2,DCM.n,1);
        DCM.TE     = 0.04;
        DCM.options.nonlinear  = 0;
        DCM.options.two_state  = 1;
        DCM.options.stochastic = 0;
        DCM.options.nograph    = 1;
        % Connectivity matrices
        DCM.a=[1 0 1 0;
            1 1 0 0;
            1 0 1 1;
            0 1 0 1];
        DCM.b = zeros(4,4,3);
        DCM.b(1,3,1) = 1;DCM.b(:,:,2) = 0;DCM.b(:,:,3) = 0;
        DCM.c =[0 0 0;
            0 0 0;
            1 0 0;
            0 1 1];
        save(fullfile(output_dir,'DCM_13.mat'),'DCM');
        % % DCM Estimation
        clear matlabbatch
        matlabbatch{1}.spm.dcm.fmri.estimate.dcmmat = {fullfile(output_dir,'DCM_13.mat')};
        spm_jobman('run',matlabbatch);
    end
    %%  Bayesian Model Comparison
    for j=1:length(sb)
        DCM_full = load('DCM_full.mat','F');
        DCM_1 = load('DCM_1.mat','F');
        DCM_2 = load('DCM_2.mat','F');
        DCM_3 = load('DCM_3.mat','F');
        DCM_4 = load('DCM_4.mat','F');
        DCM_5 = load('DCM_5.mat','F');
        DCM_6 = load('DCM_6.mat','F');
        DCM_7 = load('DCM_7.mat','F');
        DCM_8 = load('DCM_8.mat','F');
        DCM_9 = load('DCM_9.mat','F');
        DCM_10 = load('DCM_10.mat','F');
        DCM_11 = load('DCM_11.mat','F');
        DCM_12 = load('DCM_12.mat','F');
        DCM_13 = load('DCM_13.mat','F');        
        fprintf('Model evidence: %f (full) vs %f (_1)\n',DCM_full.F,DCM_1.F,DCM_2.F,DCM_3.F,DCM_4.F,DCM_5.F,DCM_6.F,DCM_7.F,DCM_8.F,DCM_9.F,DCM_10.F,DCM_11.F,DCM_12.F,DCM_13.F);
        clear matlabbatch;
        
        %BMA models comparison & review
        matlabbatch{1}.spm.dcm.bms.inference.dir = {fullfile(root,sb{j},ana_dir)};
        matlabbatch{1}.spm.dcm.bms.inference.sess_dcm{1}.dcmmat = {
            fullfile(root,sb{j},ana_dir,'DCM_full.mat')
            fullfile(root,sb{j},ana_dir,'DCM_1.mat')
            fullfile(root,sb{j},ana_dir,'DCM_2.mat')
            fullfile(root,sb{j},ana_dir,'DCM_3.mat')
            fullfile(root,sb{j},ana_dir,'DCM_4.mat')
            fullfile(root,sb{j},ana_dir,'DCM_5.mat')
            fullfile(root,sb{j},ana_dir,'DCM_6.mat')
            fullfile(root,sb{j},ana_dir,'DCM_7.mat')
            fullfile(root,sb{j},ana_dir,'DCM_8.mat')
            fullfile(root,sb{j},ana_dir,'DCM_9.mat')
            fullfile(root,sb{j},ana_dir,'DCM_10.mat')
            fullfile(root,sb{j},ana_dir,'DCM_11.mat')
            fullfile(root,sb{j},ana_dir,'DCM_12.mat')
            fullfile(root,sb{j},ana_dir,'DCM_13.mat')};
        matlabbatch{1}.spm.dcm.bms.inference.model_sp = {''};
        matlabbatch{1}.spm.dcm.bms.inference.load_f = {''};
        matlabbatch{1}.spm.dcm.bms.inference.method = 'RFX';
        matlabbatch{1}.spm.dcm.bms.inference.family_level.family_file = {''};
        matlabbatch{1}.spm.dcm.bms.inference.bma.bma_yes.bma_famwin = 'fanwin';
        matlabbatch{1}.spm.dcm.bms.inference.bma.bma_no = 0;
        matlabbatch{1}.spm.dcm.bms.inference.verify_id = 1;
        matlabbatch{2}.spm.dcm.bms.results.bmsmat(1) = cfg_dep('Model Inference: BMS.mat File', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','bmsmat'));
        spm_jobman('run',matlabbatch);
        
        %find model.post winning model per subject & its column
        [cell, col] = max(BMS.DCM.rfx.model.post);
        WinModel = [cell col];
        save ('model','WinModel');        
    end
    
end

