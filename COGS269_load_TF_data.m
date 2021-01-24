function [tf_data, p] = COGS269_load_TF_data(cfg)
% function [tf_data, p] = COGS269_load_TF_data(cfg)
%
% Generic Function that loads the Field Trip Time-Frequency structures that
% were saved previously from scripts like COGS269_TF_in_FT.m
%
% Input: cfg, a structure containing the various fields needed to load in
% the correct data and apply the correct baseline, if necessary.
%
% Outputs:
% tf_data = a cell array of FT structures containing the time-freq data.
% p = the parameter field used later on (powspctrm or itpc)
% 
% K. Backer, 1 May 2017

% Load in Data, already in FT structure, needed for running the permutation
% test, variable is called "data_out".
if cfg.an_type == 1
    load([cfg.in_dir,'Power_',cfg.vis_cond,'_n=5.mat']);
    p = 'powspctrm';
    
    % Power Data aren't baselined yet, so do that here:
    cfg_base = [];
    cfg_base.baseline = cfg.baseline;
    cfg_base.baselinetype = cfg.baselinetype; % note: changing this will affect the zlimit
    % that should be used for plotting...
    cfg_base.parameter = p;
    for x = 1:length(TFdata) % loop through each subject's freq data.
        TFdata{x} = ft_freqbaseline(cfg_base,TFdata{x});
    end % for x
elseif cfg.an_type == 2
    load([cfg.in_dir,'ITPC_',cfg.vis_cond,'_n=5_95Trials.mat']);
    p = 'itpc';
end

tf_data = TFdata;