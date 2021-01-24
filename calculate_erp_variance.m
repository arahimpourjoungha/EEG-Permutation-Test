function [out_var] = calculate_erp_variance(data, vflag)
% function [out_var] = calculate_erp_variance(data, vflag)
%
% INPUTs:
% data = cell array of the individual subject data, where each cell
% corresponds to a different experimental condition or group, and the 3D data
% matrix within each cell is arranged: Channel x Time x Subject
%
% vflag = a flag to tell the program what kind of variance measure should
% be calculated and outputted.
%
% If vflag = 1: compute within-subjects standard error (the way this is
% actually calculated depends on the number of conditions inputted. If 2
% conditions, it will subtract one from the other and compute the SEM on
% the difference wave --> each condition ends up with the same variance.
% If more than 2 conditions, it will mean center the data and then take the
% standard error of each condition's mean centered data separately.)
%
% If vflag = 2: Computes regular (across-subjects) standard error (default).  As a
% result, each condition will have a different variance associated with it,
% even if there are 2 conditions, unlike the within-subjects sem output. 
% This method will also be used if data contains only 1 cell (1 condition/group). 
%
%
% OUTPUT:
% out_var = a cell array of the variance computed, where each cell
% corresponds to a different condition or group (same # of cells as input data)
% and each cell contains a matrix(channels x times), with the
% variance values.
%
% K. Backer, November 2016.

% First check inputs to make sure all is kosher.
num_cells = length(data);
if exist('vflag','var')
    if num_cells == 1 && vflag ~= 2
        vflag = 2;
        display('Only 1 Data Cell detected. Changing Variance Flag to 2 (compute across-subjects SEM)')
    elseif vflag == 1 || vflag == 2
        display('Passed Input Variable Check. Using user-defined variance calculation method.')
    elseif vflag == 0 || vflag > 2
        error('Wrong Input value used for Variance Flag. It should be either 1 or 2.')
    end % if num_cells
else % Vflag wasn't inputted. Defaults to across-subject SEM.
    vflag = 2;
    display('No Variance Flag input. Defaulting to Across-Subjects Standard Error.')
end

% Set up the Out Variance variable:
out_var = cell(size(data));
for x = 1:length(out_var)
    out_var{x} = zeros(size(data{1},1),size(data{1},2));
end % for x

if vflag == 1 % Within-subjects SEM computation (a little more complicated)
    display('Computing Within-Subjects SEMs...')
    if num_cells == 2 % Paired t-tests assumed, so calcuate corresponding variance:
        % Both out_var cells will have the same variance values by
        % definition:
        display('2 conditions detected.  Paired t-tests assumed. Calculating SEM of difference waves.')
        diff_data = data{1} - data{2};
        for ch = 1:size(diff_data,1) % loop through each channel
            for t = 1:size(diff_data,2) % Loop through each time point
                temp_dd = squeeze(diff_data(ch,t,:));
                out_var{1}(ch,t,:) = std(temp_dd)/sqrt(length(temp_dd));
            end % for t
        end % for ch
        out_var{2} = out_var{1}; % both variances are equal.
        
    elseif num_cells > 2 % Assume ANOVA stats, use a mean-centering approach.
        % Each out_var cell will have different variance values:
        display('More than 2 conditions detected. Using mean-centering approach.')
        
        mc_data = cell(size(data)); % To Hold the Mean-Centered Data.
        for x = 1:length(mc_data)
            mc_data{x} = zeros(size(data{1}));
        end % for x
        
        for ch = 1:size(data{1},1)
            for t = 1:size(data{1},2)
                % Create a matrix to hold a subset of the data:
                % subjects x conditions.
                temp_data = zeros(size(data{1},3),num_cells);
                for d = 1:num_cells
                   temp_data(:,d) = squeeze(data{d}(ch,t,:))';                     
                end % for d
                
                % Now, find the mean of temp_data, for each subject:
                temp_mean = mean(temp_data,2); % average across columns.
                
                % Now, subtract this temp mean from temp_data:
                temp_mc_data = temp_data - repmat(temp_mean,1,num_cells);
                
                % Finally, put the mc_data back into the appropriate place
                % in mc_data:
                for d = 1:num_cells
                    mc_data{d}(ch,t,:) = temp_mc_data(:,d);
                end % for d
            end % for t
        end % for ch
        
        % Now, finally, find the SEM for each cell, ch, and time of the
        % mean-centered data:
        for d = 1:num_cells
            % Within each cell, loop through each time point (2nd dimension of
            % data), and compute the SEM for each time point:
            for ch = 1:size(mc_data{d},1)
                for t = 1:size(mc_data{d},2)
                    temp_dd = squeeze(mc_data{d}(ch,t,:));
                    out_var{d}(ch,t,:) = std(temp_dd)/sqrt(length(temp_dd));
                end % for t
            end % for ch
        end % for d
    end
    
elseif vflag == 2 % Across-subjects SEM computation (easier)
    display('Computing Across-Subjects SEMs...')
    % Simply loop through each data cell:
    for d = 1:num_cells
        % Within each cell, loop through each time point (2nd dimension of
        % data), and compute the SEM for each time point:        
        for ch = 1:size(data{d},1)
            for t = 1:size(data{d},2)
                temp_dd = squeeze(data{d}(ch,t,:));
                out_var{d}(ch,t,:) = std(temp_dd)/sqrt(length(temp_dd));                
            end % for t
        end % for ch
    end % for d    
end % if vflag