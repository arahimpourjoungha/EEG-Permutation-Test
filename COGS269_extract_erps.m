function [erp_data,erp_var,idata] = COGS269_extract_erps(cfg)
% function [erp_data,erp_var,idata] = COGS269_extract_erps(cfg)
%
% THIS LOADS IN THE RELEVANT ERPLAB DATA FILE AND THEN USES OTHER INFO IN
% CFG TO EXTRACT THE CORRECT SUBJECTS AND BINS FROM THE ERPLAB STRUCTURE.
%
% INPUT:
% CFG = A STRUCTURE, containing a variety of fields, containing all the
% information needed to load in the right data and extract the right bins.
%
% OUTPUTS:
% ERP_DATA = A CELL ARRAY OF FTSTRUCT'S, 1 PER REQUESTED SUBJECT.
% ERP_VAR = A MATRIX CONTAINING THE VARIANCE (SEE
% CALCULATE_ERP_VARIANCE.M) AT EACH TIME POINT, DEFAULT = ACROSS-SUBJECTS
% SEM.
% IDATA = A 3D MATRIX CONTAINING ALL SUBJECT'S ERPS
%
% K. BACKER, 29 OCTOBER 2019

% Setup variables to collect the ERPs:
erp_data = cell(size(cfg.binnames));
erp_var = cell(size(cfg.binnames));
for c = 1:length(erp_data)
    erp_data{c} = cell(size(cfg.reqSID));
end % for e
idata = cell(size(cfg.binnames)); % Keep individual subject data in matlab format too.

% Loop through each subject:
for x = 1:length(cfg.reqSID)
    fn = [num2str(cfg.reqSID(x)),'-ERPs.erp']; % set the filename, e.g., 1-ERPs.erp
    [ERP] = pop_loaderp('filename',fn,'filepath',cfg.erpdir);
    
    % Remove the extra channels from the ERP structure, so that the channels match that in the
    % EEGLAB structure:
    erps = ERP.bindata((1:cfg.num_chans),:,:);
    
    if x == 1 % if it's the first subject, do some checks...
        % Now, before replicating the fieldtrip structure for each condition/subject, check the FieldTrip
        % structure to make sure the times and sampling rate is correct.
        % Also, remove the variance field, since that will change, and I'm not
        % sure if it's really needed to run the permutation tests.
        ftstruct = cfg.ftstruct;
        if isfield(ftstruct,'var')
            ftstruct = rmfield(ftstruct,'var'); % remove variance field.
        end
        %if ~isequal(ftstruct.fsample,ERP.srate)
        ftstruct.fsample = ERP.srate; % verify sampling rate.
        %end
        if ~isequal(ftstruct.time,ERP.times)
            ftstruct.time = ERP.times/1000; % verify times.
            % ERP Time is in milliseconds... convert to seconds for FieldTrip.
        end
    end % if x == 1
    
    % Make sure our bins are in the correct order: 
    % Bin 1 should be Target and Bin 2 should be Standard:
    if ~strcmpi(ERP.bindescr{1}(1:6),cfg.binnames{1}(1:6)) || ...
            ~strcmpi(ERP.bindescr{2}(1:6),cfg.binnames{2}(1:6)) %check first 6 characters
        error('Binnames are in the wrong order. Fix before proceeding.')
    else % Bin indexing is correct.
        for c = 1:length(cfg.binnames) % loop through each condition.          
            tempft = ftstruct;
            tempft.avg = squeeze(erps(:,:,c)); % put in field trip structure.
            
            % now, add the new fieldtrip structure to cond_data{c}{s}:
            erp_data{c}{x} = tempft;
        end % for c
    end % if
end % for x

% Calculate the SEM of the ERP data:
if strcmpi(cfg.tt_type,'depsamplesT')
    for c = 1:length(cfg.binnames)
        % Make a 3-D matrix (idata): channels x times x subjects:
        idata{c} = zeros(size(erp_data{c}{1}.avg,1),size(erp_data{c}{1}.avg,2),length(cfg.reqSID));
        for s = 1:length(cfg.reqSID) % loop through each subject.
            idata{c}(:,:,s) = erp_data{c}{s}.avg; % add each subject's data to idata.
        end % for y
    end % for c
    % Now, calculate the standard error at each ERP time point.
    semflag = 1; % 1 for within-subjects SEM; 2 for across-subjects SEM
    [temp_var] = calculate_erp_variance(idata,semflag); % custom function.
    erp_var = temp_var;
    % Later, the Plot Code assumes data and variance is time x channels...
    % But ERPLAB structure is channels x times... need to transpose
    % the variance, as done also for the data.
    for c = 1:length(erp_var)
        if size(erp_var{c},1)== cfg.num_chans
            erp_var{c} = erp_var{c}';
        end % if
    end % for c
end