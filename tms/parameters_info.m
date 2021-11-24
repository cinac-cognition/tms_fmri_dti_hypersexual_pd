function [parameters strings] = parameters_info()
check = 1
while check
    prompt = {'Select a dataset (real / sham OR both). If both is selected, intreprete the output with caution)', ...
        'Do you want to filter outliers? If so, write how many SD (0 for none, 3.5, 3 or 2.5)', ...
        'Would you like to replace go omissions with the max RT of go trials? (yes / no)', ...
        ['Select one of the following conditions (write the condition number below): \n \n'...
        '1. - Condition 1: yields the average (i/m)SSRT using the last half of SSD trials \n'...
        '2. - Condition 2: yields the average (i/m)SSRT using all SSD trials \n'...
        '3. - Condition 3: yields the median (i/m)SSRT using the last half of SSD trials  \n'...
        '4. - Condition 4: yields the median (i/m)SSRT using all SSD trials  '...
        '5. - Condition 5: Computes all the above conditions at once. Interprete the output with caution. \n \n Your selection is: \n']}
    title='Please enter session infomation';
    boxsize=repmat([1,70],[4 1]); %set size of entry box to 70X1 (in 6x1 matrix [each dialog window]
    tmp=inputdlg([{sprintf(prompt{1})} {sprintf(prompt{2})} {sprintf(prompt{3})}  {sprintf(prompt{4})}], title, boxsize)
    
    if  (strcmp(tmp{1}, 'real')) |  (strcmp(tmp{1}, 'sham'))
        datarun = 1
    elseif  (strcmp(tmp{1}, 'both'))
        datarun = 1:2
    end
    
    %Add the info to a new structure
    parameters = struct('dataset', tmp{1}, 'datarun', datarun, 'filter', str2double(tmp{2}), 'replacement', tmp{3}, 'condition', str2double(tmp{4}));
    
    if parameters.condition == 5
        parameters.condition = 1:4
    end
    
    
    
    
    if ~(strcmp(parameters.dataset, 'real') | strcmp(parameters.dataset,'sham') |  strcmp(parameters.dataset,'both') )
        f = msgbox('Error, invalid dataset, please, re-run this code.')
        return
    end
    if (strcmp(parameters.replacement, 'yes') | strcmp(parameters.replacement,'no'))
        if (strcmp(parameters.replacement, 'yes'))
            parameters.replacement = true
        elseif (strcmp(parameters.replacement, 'no'))
            parameters.replacement = false
        end
    else
        f = msgbox('Error, invalid replacement option, please, re-run this code.')
        return
    end
    
    if ~isnumeric(parameters.filter)
        f = msgbox('Error, invalid filter value, please, re-run this code.')
        return
    end
    if ~isnumeric(parameters.filter)
        f = msgbox('Error, invalid condition value, please, re-run this code.')
        return
    end
    
    
    strings.condition.all = {' Condition 1: yields the average (i/m)SSRT using the last half of SSD trials' ...,
        ' Condition 2: yields the average (i/m)SSRT using all SSD trials' ...,
        ' Condition 3: yields the median (i/m)SSRT using the last half of SSD trials' ...,
        ' Condition 4: yields the median (i/m)SSRT using all SSD trials'}
    
    strings.filter.true = 'Trials above or below %d std were discarded from the analysis'
    strings.filter.false = 'No filter was applied and therefore no trials were discarded from the analysis'
    strings.replacement.true = 'Trials with go omission errors were replaced by the max RT of the correct go trials distribution'
    strings.replacement.false = 'Trials with go omission errors were not replaced and therefore excluded from the analysis'
    
    
    
    
    check = 0
end
end