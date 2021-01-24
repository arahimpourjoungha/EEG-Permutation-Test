function [stat,cfg] = COGS269_run_ERP_FT_stats(tt_type,neighbours,data1,data3)
% function [stat,cfg] = COGS269_run_ERP_FT_stats(tt_type,neighbours,data1,data2)
%
% THIS CODE RUNS PERMUTATION-BASED STATS IN FIELDTRIP.
% THE OPTIONS ARE SET DEPENDING MAINLY ON THE TT_TYPE (T-TEST TYPE), AS
% WELL AS THE ERP NAME, AND FOR PAIRED-SAMPLES TESTS THE CONTRAST.
%
% K. BACKER, 11 APRIL 2017

% General cfg Variables
cfg = [];
% cfg.channel: default = all
% cfg.latency: [begin end] in seconds, default = all
cfg.avgoverchan = 'no';
cfg.avgovertime = 'no';
cfg.parameter = 'avg'; % or trial
cfg.method = 'montecarlo';

% set up all the other cfg variables for method montecarlo:
cfg.numrandomization = 5000;% number of permutations
cfg.correctm = 'cluster'; % correction method.
cfg.alpha = 0.7;

cfg.randomseed = 'yes'; % default or no.
cfg.statistic = tt_type; % Type of T-Test to run:
% 'onesampleT', 'depsamplesT', 'indepsamplesT'

cfg.tail = 0; % Default.
cfg.correcttail = 'prob'; % If doing a 2-tailed t-test, prob is recommended.
% Use 'no' if doing 1-tailed test. Another way (not recommended) to correct
% is 'alpha'.

if strcmpi(cfg.correctm,'cluster')
    % since cluster was selected above, input more cfg parameters!
    cfg.clusterstatistic = 'maxsum';%'wcm'; % maxsum (default) or maxsize or wcm
    cfg.clusterthreshold = 'parametric'; % or 'nonparametric_individual' or
    % 'nonparametric_common'
    cfg.clusteralpha = 0.8; % default = 0.05
    % cfg.clustercritval --> default determined by the stat function used.
    cfg.neighbours = neighbours;
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
if strcmpi(tt_type,'depsamplesT')
    nsub = length(data1); % number of subjects, cond. 1
    nsub2 = length(data3); % number of subjects, cond. 2
    % Check to make sure they're equal:
    if nsub ~= nsub2
       error('Different number of subjects with condition 1 and condition 2 data!') 
    end
    cfg.design(1,1:2*nsub) = [ones(1,nsub) 2*ones(1,nsub)];
    cfg.design(2,1:2*nsub) = [1:nsub 1:nsub];
    cfg.ivar = 1; % Row number in design corresponding to independent variable(s)
    cfg.uvar = 2; % Row number in design corresponding to unit varaibles (i.e., subjects)
    % Unit of Observation variables... subject #s or trial #s
elseif strcmpi(tt_type,'indepsamplesT')
    nsub = length(data1); % number of subjects, group 1
    nsub2 = length(data3); % number of subjects, group 2
    % Currently, we have different numbers of subjects per group...
    % variables below have been adjusted, for this case:
    cfg.design(1,1:nsub+nsub2) = [ones(1,nsub) 2*ones(1,nsub2)];
    cfg.ivar = 1; % Row number in design corresponding to independent variable(s)
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmpi(tt_type,'depsamplesT') || strcmpi(tt_type,'indepsamplesT')
    % run the permutation tests!
    stat = ft_timelockstatistics(cfg,data1{:},data3{:});
else
    error('Inputted t-test type not recognized.')
end
