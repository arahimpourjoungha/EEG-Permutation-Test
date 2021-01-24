function [layout,cfg] = COGS269_make_FT_layout(ftstruct,num_chans,ext)
% function [layout,cfg] = COGS269_make_FTlayout_file(ftstruct,num_chans,ext)
%
% This function makes the layout file in ext (.mat or .lay) format based on the number of
% channels specified and using the electrode coordinates in ftstruct.
%
% K. Backer, 10 APRIL 2017

cfg = [];
cfg.elec = ftstruct.elec;
cfg.projection = 'polar';
cfg.output = ['Biosemi',num2str(num_chans),ext];
cfg.rotate = 90;
% not setting the next two options, leads to 2 extra "Channels", 1
% for the scale and 1 for the comment... not sure if this will
% break anything.
cfg.skipscale = 'yes';
cfg.skipcomnt = 'yes';
[layout, cfg] = ft_prepare_layout(cfg);
figure,ft_plot_lay(layout)