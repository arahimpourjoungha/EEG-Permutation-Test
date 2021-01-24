% COGS269_TF_Stats.m
%
% Main script to run stats on the time-frequency data.
%
% K. Backer, 2019 Nov 19

clear; % Clear workspace.

% Set Fixed Variables:
subjects = 1:5; % which subjects to include in the statistical analysis.
in_dir = 'C:\Users\kbacker\Documents\Teaching\COGS269_EEG\Data\TF\'; % Where your TF files are stored.
an_type = 1; % 1 for Power, 2 for ITPC

% Set Baseline Parameters for Power Analysis:
cfg = [];
cfg.baselinetype = 'db';% 'absolute', 'relative', 'relchange', 'normchange', 'db', 'vssum' or 'zscore' (default = 'absolute')
cfg.baseline = [-0.5 0]; %in seconds
cfg.parameter = 'powspctrm'; % the field in the TFdata structure to apply baselining to.

% Type of t-test: Dependent-Samples t-test to contrast Standards vs. Targets.
tt_type = 'depsamplesT';

% Addpath to fieldtrip (if not already in path):
%addpath 'C:\Users\kbacker\Documents\fieldtrip-20191025\';
ft_defaults; % Load FieldTrip defaults.

% Addpath to EEGLAB (need for copyaxis, to click on figures):
%addpath('C:\Users\kbacker\Documents\eeglab_current\eeglab2019_0\');
eeglab;
%% Load in Time-Frequency Data, Also Baseline if an_type = 1 (power).
vis_conds = {'Target' 'Standard'}; % Names of the visual conditions

% Loop through each visual condition and load in Time-Frequency (Power or
% ITPC) data:
all_tf_data = cell(size(vis_conds)); % cell array with all data.
for v = 1:length(vis_conds)
 
    cfg.an_type = an_type;
    cfg.in_dir = in_dir;
    cfg.vis_cond = vis_conds{v};
    
    cfg_base = [];
    cfg_base.baseline = cfg.baseline;
    cfg_base.baselinetype = cfg.baselinetype; % note: changing this will affect the zlimit
    % that should be used for plotting...
   
    [tf_data, p] = COGS269_load_TF_data(cfg);
    
    all_tf_data{v} = tf_data;
end 
%% Run Permutation Tests

% Set up Channel Neighborhoods, Assuming Layout File already exists!
ext = '.mat'; 
num_chans = length(tf_data{1}.label); % get number of channels.
lay_name = ['Biosemi',num2str(num_chans),ext];
[neighbours] = COGS269_setup_FT_neighbours(lay_name);

% Organize Inputs for Permutation Tests.
in_cfg = [];
in_cfg.tt_type = tt_type; % t-test type
in_cfg.an_type = an_type; % analysis type (Power or ITPC)
in_cfg.neighbours = neighbours; % neighbour map
in_cfg.p = p; % which parameter (field) in structure to run the perm tests on: powspctrm or itpc
in_cfg.data1 = all_tf_data{1};
in_cfg.data2 = all_tf_data{2};

% this is the command that actually runs the permutation tests:
% open this to view/change the settings for the permutation tests:
[stat,out_cfg] = COGS269_run_TF_FT_stats(in_cfg);

%% Plot Results

% Take the grand average of the data for each condition:
cfg2 = [];
cfg2.parameter = p;
cfg2.keepindividual = 'no';
grandtf1 = ft_freqgrandaverage(cfg2,all_tf_data{1}{:});
grandtf2 = ft_freqgrandaverage(cfg2,all_tf_data{2}{:});

% Extract the grand averaged data values from the grandtf structure, for plotting below.
granddata1 = eval(['grandtf1.', p]); % Extract the grand averaged data values from the grandtf structure.
granddata2 = eval(['grandtf2.', p]); % Extract the grand averaged data values from the grandtf structure.

% Set plotting limits for any colorbar:
if an_type == 1 % power
    if strcmpi(cfg_base.baselinetype,'relchange')
        zlimits = [-1 1]; % change this according to type of power baseline and the analysis.
    elseif strcmpi(cfg_base.baselinetype,'db')
        zlimits = [-3 3];
    elseif strcmpi(cfg_base.baselinetype,'absolute')
        zlimits = [-4 4];
    elseif strcmpi(cfg_base.baseline,'no') % No Power baseline
        zlimits = [-10 10];
    else
        zlimits = [-1 1];
    end
    
elseif an_type == 2 % iTPC
    zlimits = [-0.5 0.5]; % Note these values are always positive, but I like to set 
    % the scale from - to + so that the hot colors are all positive.
end

% Plot the Raw TF Spectrograms -- Grand Averaged across subjects:
plot_channel_TFdata({granddata1 granddata2}, num_chans, grandtf1.time, grandtf1.freq, zlimits,vis_conds);

% If cluster correction was done, plot the cluster results:
if strcmpi(stat.cfg.correctm,'cluster')
    plot_cluster_resultsTF(stat,{grandtf1 grandtf2},zlimits,p,vis_conds);
end