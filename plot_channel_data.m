function [fdr_mask] = plot_channel_data(data, num_chans, xvalues, ylimits, plot_colors, ribbon_flag, var_plot, pvals)
% function [fdr_mask] = plot_channel_data(data, num_chans, xvalues, ylimits, plot_colors, ribbon_flag, var_plot, pvals)
%
% This function will plot all the data inputted into 1 large figure.
%
% INPUTS:
% DATA = the data you want to plot. To plot multiple conditions or groups
% on the same plot, input a cell array of data matrices.
% Data matrices should be in the format [times x channels].
%
% NUM_CHANS = the number of channels (61 or 65)
% XVALUES = the times to plot along the x-axis.
% YLIMITS = the max and min values for the y-axis.
% PLOT_COLORS = a cell array of the colors you want to use in the plot;
% must match the number of data matrices inputted.
% RIBBON_FLAG, If set to 1, it will plot the variance inputted as var_plot
% in ribbons surrounding the means.
% VAR_PLOT = the variance you want to plot; should be the same size as
% data.
% PVALS = a cell array of the pvals you want to include on the plots.
%
% OUTPUT:
% FDR_MASK = if FDR correction is applied, the mask returned by EEGLAB.
%
% K. Backer, 21 April 2016

eval(['electrodechart',num2str(num_chans)]);

if ~exist('plot_colors','var')
    plot_colors = {'r' 'b' 'm' 'c' 'k' 'g'};
    %plot_colors = {'b' 'k' 'c' 'g'};
end
if ~exist('ribbon_flag','var')
    ribbon_flag = 0;
end

raster_flag = 0; % don't plot rasters here, it's done in another function.

% Check to see if ylimits is a cell array, as it is now for new CMP stats
% code... 
% If so, take the largest ranges:
if iscell(ylimits)
    y = cell2mat(ylimits);
    ylimits = [min(y) max(y)];
end

if exist('pvals','var')
    thold = 0.8; % clusters are already corrected. % 0.005; % uncorrected threshold
    fdr_flag = 0;% 0 or 1 to do FDR correction.
    fdr_thold = 0.8; % Threshold for FDR correction.
    fdr_type = 'parametric'; % ['parametric' (Def) |'nonParametric'] FDR type.
    %p_ranges = {[4.5 5] [4 4.5] [3.5 4]}; % y-axis ranges for plotting the significant time bins
    p_ranges = {[ylimits(2)-0.5 ylimits(2)] [ylimits(2)-1 ylimits(2)-0.5] [ylimits(2)-1.5 ylimits(2)-1]};
    %p_ranges2 = {[-3.5 -4] [-4 -4.5] [-4.5 -5]}; % FDR p-values
%     if ylimits(1) > min(p_ranges2{3})
%         ylimits(1) = min(p_ranges2{3});
%     end
%     if ylimits(2) < max(p_ranges{1})
%         ylimits(2) = max(p_ranges{1});
%     end
    
    % for each set of pvals entered, find the pvals that meet the uncorrected
    % threshold.  Also, if fdr_flag is set, find the pvals that are below the
    % fdr threshold.
    fdr_mask = cell(size(pvals));
    uncorr_mask = cell(size(pvals));
    for x = 1:length(pvals)
        % Make a mask for the uncorrected pvals:
        z = zeros(size(pvals{x}));
        f = find(pvals{x}<thold);
        z(f) = 1;
        uncorr_mask{x} = z;
        
        if fdr_flag == 1
            [p_fdr, p_masked] = fdr(pvals{x},fdr_thold,fdr_type);
            fdr_mask{x} = p_masked;
        end
        
        % Check to see if xvalues is longer than the masks.
        % if so, add in 0's to the beginning of each mask.
        
%         if length(xvalues) > size(fdr_mask{x},1)
%             diff_length = length(xvalues)-size(uncorr_mask{x},1);
%             zeros_to_add = zeros(diff_length,size(uncorr_mask{x},2));
%             uncorr_mask{x} = [zeros_to_add; uncorr_mask{x}];
%             if fdr_flag == 1
%                 fdr_mask{x} = [zeros_to_add; fdr_mask{x}];
%             end
%         end
    end % for x
end % if exist


figure

% Loop through each Channel:
for ch = 1:num_chans
    subplot(num_rows, num_cols, plot_idx{ch}{2});
    
    % Loop through each group or condition inputted:
    ihdls = {};
    ihdls2 = {};
    for x = 1:length(data)
        
        if ribbon_flag == 1 % Plot the variance in ribbons:
           
            ribbon1a = data{x}(:,plot_idx{ch}{3}) + var_plot{x}(:,plot_idx{ch}{3});
            ribbon1b = data{x}(:,plot_idx{ch}{3}) - var_plot{x}(:,plot_idx{ch}{3});
            ribbon1 = [ribbon1b' fliplr(ribbon1a')];
            
            ihdls2{x} = patch([xvalues fliplr(xvalues)], ribbon1,plot_colors{x});
            alpha(0.5);
            %alpha(1);
            hold on;
            
            ihdls{x} = plot(xvalues,data{x}(:,plot_idx{ch}{3}),plot_colors{x},'LineWidth',3);
            
        else
            ihdls{x} = plot(xvalues,data{x}(:,plot_idx{ch}{3}),plot_colors{x},'LineWidth',3);
            hold on;
            
        end % if
    end % for x
    
    % Plot a line going through 0.
    plot([0 0], [ylimits(1) ylimits(2)],'--k');
    plot([xvalues(1) xvalues(end)], [0 0],'--k');
    
    % Add in the tickmarks marking significant time bins.
    if exist('pvals','var')
        for p = 1:numel(uncorr_mask)
            chp = uncorr_mask{p}(:,plot_idx{ch}{3});
            for q = 1:size(chp,1)
                if chp(q) == 1
                    plot([xvalues(q) xvalues(q)],[p_ranges{p}(1)+abs(p_ranges{p}(1)*0.1) p_ranges{p}(1)],...
                        'LineWidth',3,'Color','k');%plot_colors{p});
                end % if
            end % for q
            
            if fdr_flag == 1
                chp = fdr_mask{p}(:,plot_idx{ch}{3});
                for q = 1:size(chp,1)
                    if chp(q) == 1
                        plot([xvalues(q) xvalues(q)],[p_ranges2{p}(1)+abs(p_ranges2{p}(1)*0.1) p_ranges2{p}(1)],...
                            'LineWidth',3,'Color','k');%plot_colors{p});
                    end % if
                end % for q
                
            end % if fdr_flag
        end % for p
    end % if exist
    
    for x = 1:length(ihdls)
        set(ihdls{x},'ButtonDownFcn','copyaxis');
        if ribbon_flag == 1
            set(ihdls2{x},'ButtonDownFcn','copyaxis');
        end
    end % for x
    
    set(gca,'XLim',[xvalues(1) xvalues(end)],'YLim',ylimits,'TickDir','out');
    axis([xvalues(1) xvalues(end) ylimits(1) ylimits(2)])
    title(plot_idx{ch}{1});
    %xlabel('Time (ms)')
    %ylabel('Amplitude (microV)')
    grid on
end % for ch

% Make Raster Plots showing significant time bins across channels!
if raster_flag == 1
    if exist('pvals','var')
        % Remove the pvals before time 0 before plotting:
%         tidx = find(xvalues == 0);
%         if isempty(tidx)
%             tidx = find(xvalues < 0);
%             tidx = tidx(end)+1;
%         end
        for p = 1:numel(uncorr_mask)
            if fdr_flag == 1
                % If the FDR is significant but the uncorr. thold isn't, then
                % this line below will create an error. the FDR value will be
                % plotted as 1 instead of 2. Set this situation to 1.5 to
                % visualize how often it happens.
                % pval_plot = uncorr_mask{p} + fdr_mask{p};
                pval_plot = zeros(size(fdr_mask{p}));
                for ch = 1:size(fdr_mask{p},2)
                    for y = 1:size(fdr_mask{p},1)
                        if fdr_mask{p}(y,ch)==1 && uncorr_mask{p}(y,ch)==1
                            pval_plot(y,ch) = 2;
                        elseif fdr_mask{p}(y,ch)==0 && uncorr_mask{p}(y,ch)==1
                            pval_plot(y,ch) = 1;
                        elseif fdr_mask{p}(y,ch)==1 && uncorr_mask{p}(y,ch)==0
                            pval_plot(y,ch) = 1.5;
                        end
                    end % for y
                end % for ch
            else
                pval_plot = uncorr_mask{p};
            end % if fdr_flag
            
            % Need to re-arrange the channel data in the order of front to
            % posterior for the raster plot.
            pval_plot2 = zeros(size(pval_plot));
            for ch = 1:num_chans
                pval_plot2(:,ch) = pval_plot(:,plot_idx{ch}{3});
            end % for ch
            pval_plot2 = pval_plot2(tidx:end,:); % Trim pre-stim 0's
            figure,imagesc(xvalues(tidx:end),[1:num_chans],pval_plot2',[0 2]);
            ax = gca;
            chan_labels = '';
            for ch = 1:num_chans
                chan_labels{ch} = plot_idx{ch}{1};
            end
            ax.YTick = [1:num_chans];
            ax.YTickLabel = chan_labels;
            grid on
            
        end % for p
    end % if exist
end