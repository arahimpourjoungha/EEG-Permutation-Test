function plot_channel_TFdataMask(rawdata, maskdata, num_chans, xvalues, yvalues, zlimits, pt,a)
% function plot_channel_TFdataMask(rawdata, maskdata, num_chans, xvalues, yvalues, zlimits, pt,a)
%
% This function will plot TF spectrograms held in data.
%
% INPUTS:
% RAWDATA = A 3-D MATRIX OF THE RAW DATA VALUES YOU WANT IN THE BACKGROUND.
% MASKDATA = A 3-D MATRIX OF THE RAW DATA VALUES YOU WANT IN THE
% FOREGROUND, WITH ALL OTHER VALUES MASKED OUT.
%
% NUM_CHANS = the number of channels (61 or 65)
% XVALUES = the times to plot along the x-axis.
% YVALUES = the frequencies for the y-axis.
% ZLIMITS = [begin end] of colormap limits, (could be a cell array)
%
% If zlimits is inputted as empty ([]), plots will not be plotted with the
% same limits, and the colorbar will be shown for each plot.
%
% PT (added March 13, 2017) = a plot title you want to use in addition to
% channel label, e.g., to append the group code name.
%
% K. Backer, 14 February 2017 -- developed for CMP data.
%
% Revised 11 May 2017, to plot the raw data with the significant parts of
% the image highlighted.

eval(['electrodechart',num2str(num_chans)]);

% Check to see if zlimits is a cell array.
% If so, take the largest ranges:
if iscell(zlimits)
    z = cell2mat(zlimits);
    zlimits = [min(z) max(z)];
end


figure

% Loop through each Channel:
for ch = 1:num_chans
    subplot(num_rows, num_cols, plot_idx{ch}{2});
    
    if ~isempty(zlimits)
        ihdl = imagesc(xvalues,yvalues,squeeze(maskdata(plot_idx{ch}{3},:,:)),zlimits);
        colormap jet;
        
        hold on
        ihdl2 = imagesc(xvalues,yvalues,squeeze(rawdata(plot_idx{ch}{3},:,:)),zlimits);
        colormap jet;
        
    else
        ihdl = imagesc(xvalues,yvalues,squeeze(maskdata(plot_idx{ch}{3},:,:)));
        colorbar
        colormap jet;
        
        hold on
        ihdl2 = imagesc(xvalues,yvalues,squeeze(rawdata(plot_idx{ch}{3},:,:)),zlimits);
        colormap jet;
    end
    axis xy
    hold on
    % Plot a line going through 0.
    plot([0 0], [yvalues(1) yvalues(end)],'--k');
    % And plot a line going through 2.5 (visual offset)
    plot([2.5 2.5], [yvalues(1) yvalues(end)],'--k');
    
    Amap = zeros(size(squeeze(rawdata(plot_idx{ch}{3},:,:))))+a;
    set(ihdl2,'alphadata',Amap);
    %set(ihdl2,'alphadata',A);
    set(ihdl2,'ButtonDownFcn','copyaxis');
    if iscell(pt)
        title([plot_idx{ch}{1},': ',pt{d}]);
    else
        title([plot_idx{ch}{1},': ',pt]);
    end
end % for ch