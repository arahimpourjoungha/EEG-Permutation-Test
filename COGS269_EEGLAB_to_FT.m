function [ftstruct] = COGS269_EEGLAB_to_FT(fn,erpdir,num_chans,num_chans_ex)
% function [ftstruct] = COGS269_EEGLAB_to_FT(fn,erpdir,num_chans,num_chans_ex)
%
% CONVERT EEGLAB DATASET TO FIELDTRIP FORMAT TO FACILITATE CONVERTING THE
% ERPLAB DATASET INTO FT STRUCTURE.
% 
% K. BACKER, 10 APRIL 2017

% Load in a different set, depending on the group b/c channel info differs
% across the groups!
EEG = pop_loadset('filename',fn,'filepath',erpdir);

% remove any extra, unwanted channels from EEG2, such as Mastoids, EOG:
extra_chans = [num_chans + 1:num_chans+num_chans_ex];
EEG2 = pop_select(EEG,'nochannel',extra_chans);
% Remember to also remove these channels from the ERP BinData that will be
% inserted into the FieldTrip structure!

% Convert this EEG2 structure to FieldTrip to serve as a template
% structure for the ERPLAB data:
ftstruct = eeglab2fieldtrip(EEG2,'timelockanalysis','none');
% Need to add the dimord field manually:
ftstruct.dimord = 'chan_time';

% Plot a blank topo while we're here!
figure; topoplot([],EEG2.chanlocs,'style','blank','electrodes','labelpoint','chaninfo',EEG2.chaninfo);