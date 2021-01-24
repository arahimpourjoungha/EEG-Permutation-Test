function [stat,cfg] = COGS269_run_TF_FT_stats(in_cfg)
% function [stat,cfg] = COGS269_run_TF_FT_stats(in_cfg)
%
% THIS CODE RUNS PERMUTATION-BASED STATS ON TIME-FREQUENCY REPRESENTATIONS IN FIELD TRIP.
% THE OPTIONS ARE SET DEPENDING MAINLY ON THE TT_TYPE (T-TEST TYPE),
%
% INPUT: in_cfg, a structure containing various fields with the info needed
% to 
%
% K. BACKER, 2 MAY 2017

% General cfg Variables
cfg.channel = 'all'; % default = all
cfg.latency = [0 0.8];% [begin end] in seconds, default = all
cfg.frequency = [8 13]; % Example: Look at Alpha Band only
%cfg.frequency = [3 30];%'all'; % This is what I would usually do with access 
% to a very powerful computer/cluster, but it probably will require too much RAM...
cfg.avgoverchan = 'no';
cfg.avgovertime = 'no';
cfg.avgoverfreq = 'no'; % can set this to 'yes' to average across frequencies. 
cfg.parameter = in_cfg.p; % Name of Field to do perm tests on.
cfg.method = 'montecarlo';

% set up all the other cfg variables for method montecarlo:
cfg.numrandomization = 5000; % number of permutations, 
cfg.correctm = 'cluster';
cfg.alpha = 0.05;
cfg.randomseed = 'yes'; % default or no.
cfg.statistic = in_cfg.tt_type; % Type of T-Test to run:
% 'onesampleT', 'depsamplesT', 'indepsamplesT'

ttail = 0; % Default.
tcorrecttail = 'prob'; % If doing a 2-tailed t-test, prob is recommended.
% Use 'no' if doing 1-tailed test. Another way (not recommended) to correct
% is 'alpha'.

cfg.tail = ttail;
cfg.correcttail = tcorrecttail;

if strcmpi(cfg.correctm,'cluster')
    % since cluster was selected above, input more cfg parameters!
    cfg.clusterstatistic = 'maxsum'; % maxsum (default) or maxsize or wcm
    cfg.clusterthreshold = 'parametric'; % or 'nonparametric_individual' or
    % 'nonparametric_common'
    cfg.clusteralpha = 0.05; % default = 0.05
    % cfg.clustercritval --> default determined by the stat function used.
    cfg.neighbours = in_cfg.neighbours;
    cfg.minnbchan = 2; % not in help, but in fieldtrip tutorial.
    if strcmpi(cfg.clusterstatistic,'wcm')        
        cfg.wcm_weight = 1;
        cfg.clustertail = cfg.tail; % default, or -1 or 1.  Set the same as tail.    
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DESIGN MATRIX CHANGES DEPENDING ON THE T-TEST!
% Set up the Design Matrix, also not in help, but in fieldtrip
% tutorial.
if strcmpi(in_cfg.tt_type,'depsamplesT')
    nsub = length(in_cfg.data1); % number of subjects, cond. 1
    nsub2 = length(in_cfg.data2); % number of subjects, cond. 2
    % Check to make sure they're equal:
    if nsub ~= nsub2
       error('Different number of subjects with condition 1 and condition 2 data!') 
    end
    cfg.design(1,1:2*nsub) = [ones(1,nsub) 2*ones(1,nsub)];
    cfg.design(2,1:2*nsub) = [1:nsub 1:nsub];
    cfg.ivar = 1; % Row number in design corresponding to independent variable(s)
    cfg.uvar = 2; % Row number in design corresponding to unit varaibles (i.e., subjects)
    % Unit of Observation variables... subject #s or trial #s
elseif strcmpi(in_cfg.tt_type,'indepsamplesT')
    nsub = length(data1); % number of subjects, group 1
    nsub2 = length(data2); % number of subjects, group 2
    % Currently, we have different numbers of subjects per group...
    % variables below have been adjusted, for this case:
    cfg.design(1,1:nsub+nsub2) = [ones(1,nsub) 2*ones(1,nsub2)];
    cfg.ivar = 1; % Row number in design corresponding to independent variable(s)
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmpi(in_cfg.tt_type,'depsamplesT') || strcmpi(in_cfg.tt_type,'indepsamplesT')
    % run the permutation tests!
    stat = ft_freqstatistics(cfg,in_cfg.data1{:},in_cfg.data2{:});
else
    error('Inputted t-test type not recognized.')
end
