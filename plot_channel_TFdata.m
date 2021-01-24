function plot_channel_TFdata(data, num_chans, xvalues, yvalues, zlimits, pt, zmap)
% function plot_channel_TFdata(data, num_chans, xvalues, yvalues, zlimits, pt, zmap)
%
% This function will plot TF spectrograms held in data.
%
% INPUTS:
% DATA = the data you want to plot held within a cell, e.g., {TFdata}. 
% To plot multiple conditions or groups, input a cell array of data matrices.
% This will plot each group or condition on a separate figure.
% Data matrices should be in the format [channels x freq x time].
%
% NUM_CHANS = the number of channels (61 or 65)
% XVALUES = the times to plot along the x-axis.
% YVALUES = the frequencies for the y-axis.
% ZLIMITS = [begin end] of colormap limits, (could be a cell array)
% PT = string, plot title.
% ZMAP = string, which colormap to use (default = 'jet')
% If zlimits is inputted as empty ([]), plots will not be plotted with the
% same limits, and the colorbar will be shown for each plot.
% 
% PT (added March 13, 2017) = a plot title you want to use in addition to
% channel label, e.g., to append the group code name.
%
% K. Backer, 14 February 2017 -- developed for CMP data.

eval(['electrodechart',num2str(num_chans)]);

% Check to see if zlimits is a cell array.
% If so, take the largest ranges:
if iscell(zlimits)
    z = cell2mat(zlimits);
    zlimits = [min(z) max(z)];
end


for d = 1:length(data)
    
    figure
    
    % Loop through each Channel:
    for ch = 1:num_chans
        subplot(num_rows, num_cols, plot_idx{ch}{2});
        
        if ~isempty(zlimits)
            ihdl = imagesc(xvalues,yvalues,squeeze(data{d}(plot_idx{ch}{3},:,:)),zlimits);
            if exist('zmap','var')
                colormap(zmap);
            else
                colormap jet;
            end
        else
            ihdl = imagesc(xvalues,yvalues,squeeze(data{d}(plot_idx{ch}{3},:,:)));
            colorbar
            if exist('zmap','var')
                colormap(zmap);
            else
                colormap jet;
            end
        end
        axis xy        
        hold on
        % Plot a line going through 0.
        plot([0 0], [yvalues(1) yvalues(end)],'--k');
        % And plot a line going through 2.5 (visual offset)
        plot([2.5 2.5], [yvalues(1) yvalues(end)],'--k');
        
        set(ihdl,'ButtonDownFcn','copyaxis'); 
        if iscell(pt)
            title([plot_idx{ch}{1},': ',pt{d}]);
        else
            title([plot_idx{ch}{1},': ',pt]);
        end
    end % for ch
end % for d