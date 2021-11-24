clc; clear all

[parameters strings] = parameters_info()
combined_ero = [];
combined_nero = [];

for dataruns = parameters.datarun
    
    clearvars -except parameters dataruns combined_ero combined_nero strings nstd
    
    %defines filepath and datafiles depending on stimulation condition
    
    if strcmp(parameters.dataset,'real') | strcmp(parameters.dataset,'both') & dataruns == 1
        if isunix
            DataPath = '/home/dmatam/pCloudDrive/rtms_icd/Real/';
        elseif IsWin
            DataPath = 'P:/DATA/rtms_icd/Real/';
        end
        strdisp.datarun = 'This is the output for the real group'
    elseif strcmp(parameters.dataset,'sham') | strcmp(parameters.dataset,'both') & dataruns == 2
         if isunix
            DataPath = '/home/dmatam/pCloudDrive/rtms_icd/Sham/';
        elseif IsWin
            DataPath = 'P:/DATA/rtms_icd/Sham/';
        end
        strdisp.datarun = 'This is the output for the sham group'
    end
    
    subjs = dir([DataPath '/ICD*']); % retrieves all task performance files (all subjs)
    cd(DataPath);
    
    %
    
    %condition can be either real/sham or both.
    
    
    for condition = parameters.condition
        
        strdisp.condition = strings.condition.all{condition}
        
        
        if parameters.condition == 1:4
            k = condition
        else
            k = length(condition)
        end
        
        for nsub = 1:length(subjs) %loop to initialize individual subject analysis
            
            
            C = importdata([DataPath subjs(nsub).name]); C.textdata = C.textdata(2:end,:); % importing data
            
            %here we eliminate those trials that
            %are "x" stds above or below the RT mean of go correct trials
            %as defined by our parameters.
            
            
            if parameters.filter > 0
                strdisp.filter = strings.filter.true
                
                corr_gotrials = C.data(:,3) == 1 & contains(C.textdata(:,2),"Go"); %go correct for obtaining RT and STD
                RTmean1(nsub,1) = nanmean(C.data(corr_gotrials,1)); %obtains mean for corr go trials.
                RTstd1(nsub,1) = nanstd(C.data(corr_gotrials,1)); %obtains STD for corr go trials.
                nstd = parameters.filter; % number of standard deviations above/below mean
                limitstd(nsub,1:2) = [RTmean1(nsub)+(nstd*RTstd1(nsub)) RTmean1(nsub)-(nstd*RTstd1(nsub))];
                delete_all = find(C.data(:,1) > limitstd(nsub,1) | C.data(:,1) < limitstd(nsub,2) & corr_gotrials == 1);
                
                C.data(delete_all,:) = [];
                C.textdata(delete_all,:) = [];
            else
                strdisp.filter = strings.filter.false
            end
            
            corr_gotrials = C.data(:,3) == 1 & contains(C.textdata(:,2),"Go"); %defines correct go trials
            rts = C.data(:,1); %defines vector of reaction times
            stimresponse = C.data(:,2); %defines vector of emitted responses 
            acc = C.data(:,3); %defines vector for accuracy
            staircase = C.data(:,4); %defines vector for stop-signal delay
            
            ntrials(nsub,1) = size(rts,1) %total number of trials
            
            stimtype = C.textdata(:,2); %type of stimuli  (go left/right and stop)
            stimdexterity = C.textdata(:,3); % 
            ncondition(nsub,1) = condition %number of condition 1 (check readme for clarifications). 
            ndatarun(nsub,1) = dataruns %number of dataruns (2 = both real/sham).
            
            
            %creating logicals
            
            go_trials = contains(stimtype,"Go"); %logical vecgtor for go trials
            go_ero = strcmp(stimtype,"Go"); go_nero =  strcmp(stimtype,"Go2"); %logical vectors for go erotic and nonerotic
            go_correct = go_trials == 1  & acc == 1; go_errors = contains(stimtype,"Go") & acc == 0; %logical vector correct and errors for go trials
            go_ero_correct = go_ero == 1 & acc == 1; go_ero_errors = go_ero == 1 & acc == 0; %logical vector for correct and incorrect erotic go trials
            go_nero_correct = go_nero == 1 & acc == 1; go_nero_errors = go_nero == 1 & acc == 0; %logical vector for correct and incorrect non-erotic go trials
            
            go_omissions = go_trials == 1  & rts == 0;  % go omission errors
            go_disc_left = go_trials == 1 & acc == 0 & stimresponse == 2; go_disc_right = go_trials == 1 & acc == 0 & stimresponse == 1; %go disc errors by dexterity
            go_disc_errors = go_disc_left + go_disc_right; %total number of go discrimination errors (left+right errors)
            
            nogo_trials = contains(stimtype,"SSD"); %logical vector defining nogo trials
            
            %ssd trials for erotic nogo trials
            SSD1 = strcmp(stimtype,"SSD1"); SSD2 = strcmp(stimtype,"SSD2"); SSD3 = strcmp(stimtype,"SSD3"); SSD4 = strcmp(stimtype,"SSD4");
            
            %ssd trials for non-erotic nogo trials
            SSDa = strcmp(stimtype,"SSDa"); SSDb = strcmp(stimtype,"SSDb"); SSDc = strcmp(stimtype,"SSDc"); SSDd = strcmp(stimtype,"SSDd");
            
            %sums up SSD for erotic / non-erotic trials
            nogo_ero = SSD1+SSD2+SSD3+SSD4; nogo_nero = SSDa+SSDb+SSDc+SSDd;
            
            %finds nogo corrects and errors in total and for each stimuli
            %(erotic/non-erotic) condition
            nogo_correct = nogo_trials == 1 & acc == 1; nogo_errors = nogo_trials == 1 & acc == 0;
            nogo_ero_correct = nogo_ero == 1 & acc == 1; nogo_ero_errors = nogo_ero == 1 & acc == 0;
            nogo_nero_correct = nogo_nero == 1 & acc == 1; nogo_nero_errors = nogo_nero == 1 & acc == 0;    %ssd trials for erotic nogo trials
            
            
            
            
            [vSSD1] = staircase(SSD1); [vSSD2] = staircase(SSD2); [vSSD3] = staircase(SSD3); [vSSD4] = staircase(SSD4);
            [vSSDa] = staircase(SSDa); [vSSDb] = staircase(SSDb); [vSSDc] = staircase(SSDc); [vSSDd] = staircase(SSDd);
           
            
            if condition == 1
            %calculating each SSD mean value taking only those from the second half of the vector     
                mSSD1 = nanmean(vSSD1(7:end)); mSSD2 = nanmean(vSSD2(7:end)); mSSD3 = nanmean(vSSD3(7:end)); mSSD4 = nanmean(vSSD4(7:end));
                mSSDa = nanmean(vSSDa(7:end)); mSSDb = nanmean(vSSDb(7:end)); mSSDc = nanmean(vSSDc(7:end)); mSSDd = nanmean(vSSDd(7:end));
            elseif condition == 2
                %calculating each SSD mean value taking entire SSD vector         
                mSSD1 = nanmean(vSSD1(1:end)); mSSD2 = nanmean(vSSD2(1:end)); mSSD3 = nanmean(vSSD3(1:end)); mSSD4 = nanmean(vSSD4(1:end));
                mSSDa = nanmean(vSSDa(1:end)); mSSDb = nanmean(vSSDb(1:end)); mSSDc = nanmean(vSSDc(1:end)); mSSDd = nanmean(vSSDd(1:end));
            elseif condition == 3
            %calculating each SSD median value taking only those from the second half of the vector
                mSSD1 = nanmedian(vSSD1(7:end)); mSSD2 = nanmedian(vSSD2(7:end)); mSSD3 = nanmedian(vSSD3(7:end)); mSSD4 = nanmedian(vSSD4(7:end));
                mSSDa = nanmedian(vSSDa(7:end)); mSSDb = nanmedian(vSSDb(7:end)); mSSDc = nanmedian(vSSDc(7:end)); mSSDd = nanmedian(vSSDd(7:end));
            elseif condition == 4
                %calculating each SSD median value taking entire SSD vector         
                mSSD1 = nanmedian(vSSD1(1:end)); mSSD2 = nanmedian(vSSD2(1:end)); mSSD3 = nanmedian(vSSD3(1:end)); mSSD4 = nanmedian(vSSD4(1:end));
                mSSDa = nanmedian(vSSDa(1:end)); mSSDb = nanmedian(vSSDb(1:end)); mSSDc = nanmedian(vSSDc(1:end)); mSSDd = nanmedian(vSSDd(1:end));
            end
            
            
            %probabily of inhibition (p|inhib) for erotic trials
            nSSD_correct_ero = sum(nogo_ero_correct);
            nSSD_failed_ero = sum(nogo_ero_errors);
            p_ero = nSSD_correct_ero / (nSSD_correct_ero + nSSD_failed_ero) * 100;
            
            %probabily of inhibition (p|inhib) for non-erotic trials
            nSSD_correct_nero = sum(nogo_nero_correct);
            nSSD_failed_nero = sum(nogo_nero_errors);
            p_nero = nSSD_correct_nero / (nSSD_correct_nero + nSSD_failed_nero) * 100;
            
            %%Calculations of iSSRT
            
            %replace go omissions with max RT value (this option is defined
            %in parameters)
            
            if parameters.replacement == true;
                strdisp.replacement = strings.replacement.true
                if ~isempty(go_omissions);
                    rts(go_omissions,1) = max(rts(find(go_correct)));
                end
                rank_ero = sort(rts(go_ero)); % this includes go_omission errors at max rt values and go discrimination errors. %sorts the distribution of rts for go erotic
                rank_nero = sort(rts(go_nero)); %sorts the distribution of rts for go non-erotic
                Ngo_ero = sum(go_ero); Ngo_nero = sum(go_nero);
                trial_i = (Ngo_ero*(p_ero+0.000001)/100); trial_i_nero = (Ngo_nero*(p_nero+0.000001)/100);
                trial_i2 = round(trial_i); trial_i2_nero = round(trial_i_nero); %deletes decimals to help search row later
            elseif parameters.replacement == false;
                strdisp.replacement = strings.replacement.false
                rank_ero = sort(rts(go_ero_correct));
                rank_nero = sort(rts(go_nero_correct));%this excludes go omission and discrimination errors from the rank.
                Ngo_ero = sum(go_ero_correct); Ngo_nero = sum(go_nero_correct); 
                trial_i = (Ngo_ero*(p_ero+0.000001)/100); trial_i_nero = (Ngo_nero*(p_nero+0.000001)/100);
                trial_i2 = round(trial_i); trial_i2_nero = round(trial_i_nero); %deletes decimals to help search row later
            end
            
            if (trial_i2 ~= 0);
                newiGo(nsub,1) = rank_ero(trial_i2,:);
                newiGo_nero(nsub,1) = rank_nero(trial_i2_nero,:); %gives back the cell value of row number from trial_i
            else
                newiGo(nsub,1) = 0;
                newiGo_nero(nsub,1) = 0;
            end
            
            %%% taks measurements
            
            %RT erotic / non-erotic with means/medians
            
            RT_ero(nsub,1) = nanmean(rts(go_ero_correct)); RT_nero(nsub,1) = nanmean(rts(go_nero_correct));
            RTmedian_ero(nsub,1) = nanmedian(rts(go_ero_correct)); RTmedian_nero(nsub,1) = nanmedian(rts(go_nero_correct));
            
            % this calculates de RT of those nogo trials that were (incorrectly)
            % responded
            
            stopRTero(nsub,1) = nanmean(rts(nogo_ero_errors));
            stopRTnero(nsub,1) = nanmean(rts(nogo_nero_errors));
            
            %retrieves all SSD staircases
            
            All_SSD_ero(nsub,:) = [mSSD1 mSSD2 mSSD3 mSSD4]; All_SSD_nero(nsub,:) = [mSSDa mSSDb mSSDc mSSDd]; %4 staircases
            
            %Average staircase SSD for erotic /n non-erotic
            SSD_ero(nsub,:) = mean(All_SSD_ero(nsub,:),2); SSD_nero(nsub,:) = mean(All_SSD_nero(nsub,:),2); % Average SSD across 4 staircases
            
            
            %calculates SSRT using classical method (assuming ~50% of prob
            %of inhibition) or iSSRT (integrative method preferred for
            %strong deviations from ~50% prob of inh) or mSSRT (sing
            %classical method with median instead of mean).
           
            SSRT(nsub,1) = RT_ero(nsub) - SSD_ero(nsub); SSRT_nero(nsub,1) = RT_nero(nsub) - SSD_nero(nsub);
            iSSRT(nsub,1) = newiGo(nsub) - SSD_ero(nsub); iSSRT_nero(nsub,1) = newiGo_nero(nsub) - SSD_nero(nsub);
            mSSRT(nsub,1) = RTmedian_ero(nsub) - SSD_ero(nsub); mSSRT_nero(nsub,1) = RTmedian_nero(nsub) - SSD_nero(nsub);
            
            
            %outputs the number of errors for each subject 
            
            ngo_disc_left_ero(nsub,1) = length(find(go_disc_left & go_ero)); ngo_disc_left_nero(nsub,1) = length(find(go_disc_left & go_nero));
            ngo_disc_right_ero(nsub,1) = length(find(go_disc_right & go_ero)); ngo_disc_right_nero(nsub,1) = length(find(go_disc_right & go_nero));
            ngo_disc_errors(nsub,1) = sum(go_disc_errors);
            ngo_errors_ero(nsub,1) = sum(go_ero_errors); ngo_errors_nero(nsub,1)= sum(go_nero_errors);
            
            %outputs the individual prob of inhibition for each subject for
            %each type of stimuli
            p_inh_ero(nsub,1) =  p_ero; p_inh_nero(nsub,1) = p_nero;
            
        end
        
        %creates final output
        Output_final = [ndatarun ncondition  RT_ero stopRTero SSD_ero All_SSD_ero SSRT iSSRT mSSRT p_inh_ero ngo_errors_ero ngo_disc_right_ero ngo_disc_left_ero...
            RT_nero stopRTnero SSD_nero All_SSD_nero SSRT_nero iSSRT_nero mSSRT_nero p_inh_nero ngo_errors_nero ngo_disc_right_nero ngo_disc_left_nero ntrials];
        
        %indicates column label
        cHeader = {'dataset' 'condition' 'RT_ero' 'stopRTero' 'SSD_ero' 'SSD1' 'SSD2' 'SSD3' 'SSD4' 'SSRT' 'iSSRT' 'mSSRT' 'p_inh_ero' 'ngo_errors_ero' 'ngo_disc_right_ero' 'ngo_disc_left_ero'...
            'RT_nero' 'stopRTnero' 'SSD_nero' 'SSDa' 'SSDb' 'SSDc' 'SSDd' 'SSRT_nero' 'iSSRT_nero' 'mSSRT_nero' 'p_inh_nero' 'ngo_errors_nero' 'ngo_disc_right_nero' 'ngo_disc_left_nero' 'ntrials'}; %dummy header
        commaHeader = [cHeader;repmat({','},1,numel(cHeader))]; %insert commaas
        commaHeader = commaHeader(:)';
        textHeader = cell2mat(commaHeader);
        %writes values to csv
        fid = fopen(['output_' 'datarun' num2str(dataruns) 'con_' num2str(condition) '.csv'],'w');
        fprintf(fid,'%s\n',strdisp.condition)
        fprintf(fid,'%s\n',strdisp.datarun)
        fprintf(fid,'%s\n',strdisp.replacement)
        fprintf(fid,'%s\n',strdisp.filter, parameters.filter)
        fprintf(fid,'%s\n',textHeader)
        fclose(fid);
        dlmwrite(['output_' 'datarun' num2str(dataruns) 'con_' num2str(condition) '.csv'],Output_final,'-append');
        
        %write data to end of file        
        parameters_analysis = parameters
        
        save(['Output_' 'datarun_' num2str(dataruns) '_con_' num2str(condition)]);
        save('parameter_analysis');
        
        
        %saves the different forms of calculating SSRT for each stimuli (erotic/non erotic) separately in a mat fil
        
        index_position = k-1;
        combined_ero{dataruns}(:,index_position+k:index_position+k+1) = [SSRT iSSRT]
        combined_nero{dataruns}(:,index_position+k:index_position+k+1) = [SSRT_nero iSSRT_nero]
        
        save('combined_ero')
        save('combined_nero')
        
        
    end
    
    
end

%Display the parameters of the analysis

if dataruns == 1:2
    
    disp('This is the ttest comparing erotic SSRT and iSSRT for real and sham groups:');
    disp(num2str([ttest2(combined_ero{1}, combined_ero{2})]));
    disp('zero means there are no significant differences between groups')
    
end