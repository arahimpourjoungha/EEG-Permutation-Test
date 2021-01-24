% COGS269_ERP_stats.m
%
% (Relatively) simple code to load in .ERP files, convert them to Field
% Trip format and run permutation stats on the ERPs.
%
% K. Backer, 27 October 2019

clear; % Clear workspace.

% Set Fixed Variables:
subjects = 1:9; % which subjects to include in the statistical analysis.
erpdir = 'C:\Users\Ali_Rahi\Desktop\EEG Course\project\Data_EEG_fingertapping\Response_Locked_ERPs\STATS'; % Where your .erp files are stored.

tt_type = 'depsamplesT'; % Type of t-test: Dependent-Samples t-test.
num_chans = 32; % include the 1st 64 (scalp) channels.
num_chans_ex = 0; % Exclude the 4 EOG and 2 mastoid channels from stats.

% Set path to EEGLAB (if not already in path):
addpath('C:\Users\Ali_Rahi\Desktop\eeglab\eeglab2019_0\');
eeglab;

% Addpath to fieldtrip (if not already in path):
addpath 'C:\Users\Ali_Rahi\Desktop\EEG Course\project\Data_EEG_fingertapping\fieldtrip-20191206\';
ft_defaults; % Load FieldTrip defaults.

% LOAD EEGLAB DATASET AND CONVERT IT TO FT STRUCTURE:
fn = '1_erp_preproc.set'; % put this .set and .fdt in the erpdir specified above.
% this eeglab dataset will be used to create the FieldTrip data structure.
[ftstruct] = COGS269_EEGLAB_to_FT(fn,erpdir,num_chans,num_chans_ex); % Custom code, call to FT

% MAKE THE LAYOUT FILE IF IT DOESN'T EXIST:
% Need it for cluster correction and topoplot.
ext = '.mat'; % save as .mat or .lay format
if ~exist(['Biosemi',num2str(num_chans),ext],'file') % Check to see if the layout file exists.
    [layout,cfg] = COGS269_make_FT_layout(ftstruct,num_chans,ext); % Custom code, call to FT to make layout
    lay_name = cfg.output;
else
    lay_name = ['Biosemi',num2str(num_chans),ext]; % It exists, no need to make it.
end % if ~exist

% now that the data is in FieldTrip structure and the channel layout has been made,
% set up the channel neighborhoods for cluster-corrected permutation tests:
[neighbours] = COGS269_setup_FT_neighbours(lay_name); % Custom code, call to fieldtrip

% Get ready to load in the ERPs:
% Which ERPLAB binnames do we want to analyze?
binnames = {'Synchronization_Pacing_Response' 'Syncopation_Pacing_Response' 'Synchronization_Continuation_Response' 'Syncopation_Continuation_Response'}; % just the first 2.
% Now, get ready to load in the .erp files for each subject and extract the
% proper bins; put all input variables into a structure:
cfg = [];
cfg.reqSID = subjects; % list of subjects
cfg.erpdir = erpdir; % the directory with the .erp files
cfg.num_chans = num_chans; % the number of channels
cfg.ftstruct = ftstruct; % Sample fieldtrip structure, made above, see COGS269_EEGLAB_to_FT.
cfg.binnames = binnames; % binnames, set above.
cfg.tt_type = tt_type; % T-Test Type.
% Now, call KB's custom code to extract ERPs from the .erp files:
[erp_data,erp_var,idata] = COGS269_extract_erps(cfg);



%% Now, Run Stats using Field Trip (finally!)
% Run stats. This function contains the call to FT to run the stats,
% and it sets up the cfg variables for FT, based on the tt_type and
% ERPn!
c1 = 1;
c2 = 2;
data1 = erp_data{c1}; % Condition 1 Data (Target)
data2 = erp_data{c2}; % Condition 2 Data (Standard)
idata1 = idata{c1};
idata2 = idata{c2};
%stat = COGS269_run_ERP_FT_stats(tt_type,neighbours,erp_data{1},erp_data{3});
stat = COGS269_run_ERP_FT_stats(tt_type,neighbours,data1,data2);

% Use the commented code below if you want to save the stat file and just load it in,
% rather than re-running the analysis. (And comment out the line above,
% starting with "stat = ...")
outfn = [erpdir,'STAT_PairedT_Targets_vs_Standards.mat']; % Filename for saving the Stats output from Field Trip
if ~exist(outfn,'file')% If the stats file doesn't exist, run stats:
   stat = COGS269_run_ERP_FT_stats(tt_type,neighbours,erp_data{1},erp_data{2});
    save(outfn,'stat');
 else % If the stats have been previously run, load that file.
     load(outfn);
 end

%% Calculate Grand Average Data in Field Trip and then Plot Grand Averaged ERPs and Results:
%data1 = erp_data{1}; % Condition 1 Data (Target)
%data2 = erp_data{2}; % Condition 2 Data (Standard)
%data3 = erp_data{3};
%data4 = erp_data{4};
% Plot the results in FieldTrip:
% First calculate the grand average for each condition:
cfg2 = [];
cfg2.channel = 'all';
cfg2.latency   = 'all';
cfg2.parameter = 'avg';
cfg2.keepindividual = 'no';
cfg2.method = 'across';
GA1 = ft_timelockgrandaverage(cfg2,data1{:}); % Grand Average 1
GA2 = ft_timelockgrandaverage(cfg2,data2{:}); % Grand Average 2
%GA3 = ft_timelockgrandaverage(cfg2,data3{:}); % Grand Average 2
%GA4 = ft_timelockgrandaverage(cfg2,data4{:}); % Grand Average 2

% This is my own plotting function.
ylimits = [-10 10]; % y-axis limits for ERP plots
gplot_colors = {'b' 'c'}; % Colors you want to use for each condition.
plot_channel_data({GA1.avg' GA2.avg'}, num_chans, GA1.time,...
    ylimits, gplot_colors, 1, erp_var, {stat.prob'});

% Plot Results of Cluster-Corrected Permutation Tests from Field Trip:
if strcmpi(stat.cfg.correctm,'cluster')
    zlimits = [-5 5]; % colorbar limits for topoplots
    %plot_cluster_results(stat,{GA1 GA2},binnames,'VEP',...
    %    {idata{1} idata{2}},ylimits,zlimits);
    plot_cluster_results(stat,{GA1 GA2},binnames,'VEP',...
        {idata1 idata2},ylimits,zlimits);
else
    warning('I have not coded results plotting for other forms of correction yet.')
end