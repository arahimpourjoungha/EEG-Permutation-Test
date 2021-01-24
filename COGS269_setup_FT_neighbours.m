function [neighbours] = COGS269_setup_FT_neighbours(lay_name)
% function [neighbours] = COGS269_setup_FT_neighbours(lay_name)
%
% This function setups the channel neighborhoods for cluster-corrected
% permutation tests.
%
% K. BACKER, 10 APRIL 2017

cfg = [];
cfg.method = 'triangulation'; 
cfg.layout = lay_name;
neighbours = ft_prepare_neighbours(cfg); % calls fieldtrip function
% Examine the neighbourhoods:
ft_neighbourplot(cfg); % Always inspect.  
% You can add or remove neighbours manually, for example (with a 32-channel
% cap), add FP1 and AF4 as neighbours:
% FP1 = Channel 1
%neighbours(1).neighblabel = {'AF3'; 'F7'; 'Fp2'; 'Fz'; 'AF4'};
% AF4 = Channel 29
% neighbours(29).neighblabel = {'AF3'; 'FC2'; 'F4'; 'F8'; 'Fp2'; 'Fz'; 'Fp1'};