function plot_cluster_results(stat,grand_erp,binnames,ERPn,idata,ylimits,zlimits)

% function plot_cluster_results(stat,grand_erp,binnames,ERPn,idata,ylimits,zlimits)
%
% INPUTS:
% stat = the stats output structure from FieldTrip.
% grand_erp = cell array of ERP data (in FT structure). The number of cells will determine
%   how many topographies and time waveforms will be plotted.
%   Input 1 cell for plotting results of 1 condition (1 sample t-test).
%   Input 2 cells for plotting results of 2 conditions (paired-samples or
%   independent samples t-test)
%   Input 6 cells for plotting results of 2x2 interaction t-test... which is a
%   either a paired-samples or independent-samples t-test done on the
%   difference waves.  The first two cells should be the difference wave
%   data, the last four cells should be the data from each of the 4
%   condition &/or group data used to make the difference waves.
%   Code doesn't currently support any other number of cells.
% binnames = cell array of strings containing the bin/condition names, for
%   labeling plots.
% ERPn = string of the ERP response being plotted, for labeling plots.
% idata = cell array of 3D matrices, containing the individual ERPs [ch x time x subs],
%   for calculating and plotting the variance as the height of the ribbons for
%   the time waveform plots
% ylimits = the ylimits for plotting time waveforms
% zlimits = the colorbar limits for plotting ERP topographies.
%
% Uses the stat output from fieldtrip to plot several things:
% 1) Raster Plot of the significant clusters.
% 2) For each cluster, topography of the T-statistic for each cluster
%   (averaged over time)
% 3) For each cluster, topography of the ERP amplitude for each cluster
%   (averaged over time); multiple topos will be plotted based on the number
%   of cells in grand_erp.
% 4) For each cluster, average time waveforms of the ERP, collapsed across
%   all channels (or top 15 channels) in the cluster.  The number of lines/plots will be determined
%   based on the number of cells in grand_erp.
%
% K. Backer, 13 December 2016.
% Revised 8 February 2017 -- Slightly different way of picking top sig
% channels to include in the cluster waveforms... now taking into account
% the direction of the effect.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PART 1. RASTER PLOT OF SIGNIFICANT CLUSTERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~isempty(strfind(ERPn,'ASSR')) || ~isempty(strfind(ERPn,'SSVEP'))
    grand_erp{1}.time = grand_erp{1}.time * 1000;
end

% custom colormap:
% light grey = n.s.
% positive clusters = warm colors
% negative clusters = cool colors
% enough colors here for 6 clusters: 3+ and 3-... if there are more
% add more colors to the custom colormap.
map = [0,1,0.7; 0,1,0; 0,0.5,0; 0,0,1; 0,0.5,1; 0,1,1; 0,0.5,0.5; 0.7,0.7,0.7; ...
    0.5,0.5,0; 1,0.5,0; 1,0,0; 0.5,0,0.5; 0.75,0,0.75;1,0,1; 1,1,0];

% Colors for Lines (for Time Waveforms):
if length(grand_erp) ~= 1
    %plot_colors = {'r' 'b' 'm' 'c' 'k' 'g'}; % Convert to map because
    %otherwise, have an error because in one case, it's referring a matrix,
    %and in another a cell.
    plot_colors = [1,0,0; 0,0,1; 1,0,1; 0,1,1; 0,0,0; 0,1,0];
else
    plot_colors = map;
end

if iscell(ylimits) 
    y = cell2mat(ylimits);
    ylimits = [min(y) max(y)];
end

if iscell(zlimits)
    z = cell2mat(zlimits);
    zlimits = [min(z) max(z)];
end

% Add some code below so that the cluster colors somewhat match the time
% waveform colors.
color_idx = [];
pos_col_idx = [9:15];
neg_col_idx = fliplr([1:7]);
zero_col_idx = pos_col_idx(1)-1;
rplot = zeros(size(stat.mask));
ppv = 1;
npv = -1;
%thold = stat.cfg.alpha; % set threshold for plotting results from parameter in stats.
thold = 0.05; % Or set a manaul threshold...
sc_idxs = {}; % Save the indices for the significant clusters for later plotting.
sc_cnt = 1;
ftstruct = grand_erp{1};
% find the p-values of the significant postive and/or negative
% clusters:
if isfield(stat,'posclusters')
    for x = 1:length(stat.posclusters) % positive clusters.
        if stat.posclusters(x).prob < thold
            % Use stat.posclusterslabelmat to get the cluster indices.
            f = find(stat.posclusterslabelmat==x);
            rplot(f) = ppv;
            color_idx(sc_cnt) = pos_col_idx(ppv);
            ppv = ppv + 1;
            sc_idxs{sc_cnt} = f;
            sc_cnt = sc_cnt + 1;
        end % if
    end % for x
end % if isfield

if isfield(stat,'negclusters')
    for x = 1:length(stat.negclusters) % negative clusters.
        if stat.negclusters(x).prob < thold
            % Use stat.negclusterslabelmat to get the cluster indices.
            f = find(stat.negclusterslabelmat==x);
            rplot(f) = npv;
            color_idx(sc_cnt) = neg_col_idx(abs(npv));
            npv = npv - 1;
            sc_idxs{sc_cnt} = f;
            sc_cnt = sc_cnt + 1;
        end % if
    end % for x
end % if isfield

% Now, re-order the channels (rows) according to the electrode
% chart and plot:
num_chans = size(stat.prob,1);
eval(['electrodechart',num2str(num_chans)]); % get channel info.
rplot2 = zeros(size(rplot));
for ch = 1:num_chans
    rplot2(ch,:) = rplot(plot_idx{ch}{3},:);
end % for ch
figure,imagesc(ftstruct.time,[1:num_chans],rplot2,[-1*((size(map,1)-1)/2) (size(map,1)-1)/2]);
ax = gca;
chan_labels = '';
for ch = 1:num_chans
    chan_labels{ch} = plot_idx{ch}{1};
end
ax.YTick = [1:num_chans];
ax.YTickLabel = chan_labels;
grid on
%colorbar
colormap(map)
if length(binnames) == 1
    title(['Significant Difference(s): ',binnames{1},' bin, for response ',ERPn]);
elseif length(binnames) == 4 % Interaction
    title(['Significant Difference(s): Interaction, for response ',ERPn]);
else
    title(['Significant Difference(s): ',binnames{1},' vs ',binnames{2},' for response ',ERPn]);
end
if ~isempty(strfind(ERPn,'ASSR')) || ~isempty(strfind(ERPn,'SSVEP'))
    xlabel('Frequency (Hz)')
else
    xlabel('Time (s)');
    hold on
    plot([0 0],[0 num_chans],':k','LineWidth',2);
end


% At this point, use the significant cluster indices to get the time and
% channels involved in each cluster.
% because the indices are not in 2D format, make a matrix of channel
% numbers the same size as stat.prob and a matrix of times the same size as
% stat.prob, using repmat.
num_times = size(ftstruct.avg,2);
chan_mat = repmat([1:num_chans]',1,num_times);
time_mat = repmat([1:num_times],num_chans,1);
sig_chans = cell(size(sc_idxs));
sig_times = cell(size(sc_idxs));
for x = 1:length(sc_idxs)
    sig_chans{x} = unique(chan_mat(sc_idxs{x}));
    sig_times{x} = unique(time_mat(sc_idxs{x}));
end % for x
% Example, for LLR, there are 2 clusters that are significant in every
% channel... Maybe add something here to take the top 15 channels based on
% the T-statistics, if all channels, or more than 15, are selected.
% Use only these top channels for plotting the ERP time waveforms for the
% clusters.
top_sig_chans = cell(size(sig_chans));
%top_chans = 16;
top_chans = floor(num_chans/2);
tstats = cell(size(sc_idxs));
for x = 1:length(top_sig_chans)
    if length(sig_chans{x}) > top_chans
        tstats{x} = stat.stat(sig_chans{x},sig_times{x});
        mean_tstat = mean(tstats{x},2); % average across time.
        %abs_mean_tstat = abs(mean_tstat); % Want the biggest values, even if negative, for sorting.
        %[s sidx] = sort(abs_mean_tstat,1,'descend');
        if color_idx(x) > zero_col_idx % Positive cluster
            [s sidx] = sort(mean_tstat,1,'descend');
        elseif color_idx(x) < zero_col_idx % Negative cluster
            [s sidx] = sort(mean_tstat,1,'ascend');
        end
        top_sig_chans{x} = sig_chans{x}(sidx(1:top_chans));
    else
        top_sig_chans{x} = sig_chans{x};
    end
end % for x

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PART 2. T-STAT TOPOS FOR EACH SIGNIFICANT CLUSTER (AVG ACROSS TIME)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
for x = 1:length(sc_idxs) % loop thru each cluster.
    subplot(length(sc_idxs),1,x)
    
    %pm = zeros(1,num_chans);
    %pm(sig_chans{x}) = 1;
    % Can't get topoplot to plot markers for only some channels, while
    % still plotting the data from those channels.
    %topoplot(mean_tstat,EEG.chanlocs,'maplimits',[-6 6]);
    
    % Try Fieldtrip instead:
    % Put the data that needs to be plotted into a temp FT structure:
    tftstruct = ftstruct;
    tftstruct.stat = stat.stat;
    % Config options:
    cfg = [];
    cfg.parameter = 'stat';
%     if strcmpi('VEP',ERPn) && x == 1 % Break the frontal positivity up to highlight the Visual P2:
%         cfg.xlim = [tftstruct.time(81) tftstruct.time(104)]; % times to average across.
%     else
%         cfg.xlim = [tftstruct.time(sig_times{x}(1)) tftstruct.time(sig_times{x}(end))]; % times to average across..
%     end
    cfg.xlim = [tftstruct.time(sig_times{x}(1)) tftstruct.time(sig_times{x}(end))]; % times to average across.
    cfg.zlim = [-8 8]; % controls colorbar limits
    cfg.marker = 'on';
    cfg.markersymbol = 'o';
    cfg.markersize = 2;
    cfg.highlight = 'on';
    cfg.highlightchannel = ft_channelselection(sig_chans{x},tftstruct.label);
    cfg.highlightsymbol = '*';
    cfg.highlightsize = 4;
    cfg.layout = ['Biosemi',num2str(num_chans),'.mat'];
    cfg.comment = 'xlim';
    cfg.commentpos = 'title';
    cfg.colorbar = 'no';
    ft_topoplotER(cfg,tftstruct);
    colormap jet;
end % for x

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PART 3. ERP AMPLITUDE TOPOS FOR EACH SIGNIFICANT CLUSTER (AVG ACROSS TIME)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set up the general config options for plotting topos of the actual data:
cfg = [];
cfg.parameter = 'avg';
cfg.marker = 'on';
cfg.markersymbol = 'o';
cfg.markersize = 2;
cfg.highlight = 'on';
cfg.highlightsymbol = '*';
cfg.highlightsize = 4;
cfg.layout = ['Biosemi',num2str(num_chans),'.mat'];
cfg.comment = 'xlim';
cfg.commentpos = 'title';
cfg.colorbar = 'no';
cfg.zlim = zlimits;

% Now, check size of grand_erp to see how many subplots to make:
if length(grand_erp) == 1 || length(grand_erp) == 2 % Only 1 condition (1-sample t-test)
    % OR 2 conditions found (paired-samples or independent samples t-test)
    sub_cnt = 1;
    figure
    for x = 1:length(sc_idxs) % loop through each cluster
        for y = 1:length(grand_erp) % loop through each cell of grand_erp (i.e., each condition)
            subplot(length(sc_idxs),length(grand_erp),sub_cnt);
%             if strcmpi('VEP',ERPn) && x == 1 % Break the frontal positivity up to highlight the Visual P2:
%                 cfg.xlim = [grand_erp{y}.time(81) grand_erp{y}.time(104)];
%             else
            cfg.xlim = [grand_erp{y}.time(sig_times{x}(1)) grand_erp{y}.time(sig_times{x}(end))]; % times to average across.
%             end
            cfg.highlightchannel = ft_channelselection(sig_chans{x},grand_erp{y}.label);
            ft_topoplotER(cfg,grand_erp{y});
            colormap jet;
            sub_cnt = sub_cnt + 1;
        end
    end
    
elseif length(grand_erp) == 6 % 6 conditions found
    % make 2 figures, 1 with difference wave topos and 1 with the original erp data from each condition
    sub_cnt = 1;
    figure % Difference wave Topo's
    for x = 1:length(sc_idxs) % loop through each cluster
        for y = 1:2 % loop through each cell of grand_erp (i.e., each condition)
            subplot(length(sc_idxs),2,sub_cnt);
            cfg.xlim = [grand_erp{y}.time(sig_times{x}(1)) grand_erp{y}.time(sig_times{x}(end))]; % times to average across.
            cfg.highlightchannel = ft_channelselection(sig_chans{x},grand_erp{y}.label);
            ft_topoplotER(cfg,grand_erp{y});
            colormap jet;
            sub_cnt = sub_cnt + 1;
        end
    end
    
    
    for x = 1:length(sc_idxs) % loop through each cluster
        sub_cnt = 1;
        figure % Condition Topo's
        for y = 3:6 % loop through each cell of grand_erp (i.e., each condition)
            subplot(2,2,sub_cnt);
            cfg.xlim = [grand_erp{y}.time(sig_times{x}(1)) grand_erp{y}.time(sig_times{x}(end))]; % times to average across.
            cfg.highlightchannel = ft_channelselection(sig_chans{x},grand_erp{y}.label);
            ft_topoplotER(cfg,grand_erp{y});
            colormap jet;
            sub_cnt = sub_cnt + 1;
        end
    end
    
    
else
    error('The number of cells in grand_erp is not currently supported.  See help file.')
end % if

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PART 4. ERP TIME WAVEFORMS FOR EACH SIGNIFICANT CLUSTER (AVG ACROSS CHAN)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xvalues = grand_erp{1}.time;
if length(grand_erp) == 1
    pc_idx = color_idx;
else
    pc_idx = [1:length(grand_erp)];
end

if length(grand_erp) == 1 || length(grand_erp) == 2 % Only 1 condition (1-sample t-test)
    % OR 2 conditions found (paired-samples or independent samples t-test)
    num_it = length(grand_erp);
elseif length(grand_erp) == 6 % Interaction
    num_it = 2;
end
for x = 1:length(top_sig_chans)
    figure
    ihdl1 = {};
    ihdl2 = {};
    %for y = 1:length(grand_erp)
    for y = 1:num_it
        % Pull out the relevant channels from grand_erp for this
        % condition and average (for the mean on the graphs)
        temp_ge = grand_erp{y}.avg(top_sig_chans{x},:);
        mean_ge = mean(temp_ge,1);
        
        % Also pull out the relevant channels from idata for this
        % condition and put into calculate_erp_variance to get the
        % variance to plot
        temp_idata = idata{y}(top_sig_chans{x},:,:);
        temp_iavg = mean(temp_idata,1); % average across the selected channels.
        [var_plot] = calculate_erp_variance({temp_iavg}, 2);
        % flag set to 2... Does across-subjects SEM for variance,
        % calculated on the individual subjects' ERP averaged across
        % the selected channels.
        
        % Now, plot the ribbon for this cluster:
        ribbon1a = mean_ge + var_plot{1};
        ribbon1b = mean_ge - var_plot{1};
        ribbon1 = [ribbon1b fliplr(ribbon1a)];
        
        if length(grand_erp) > 1
            ihdl2{y} = patch([xvalues fliplr(xvalues)], ribbon1,plot_colors(pc_idx(y),:));
            alpha(0.5);
            hold on;
            % This plots a Line for the Mean
            ihdl1{y} = plot(xvalues,mean_ge,'Color',plot_colors(pc_idx(y),:),'LineWidth',3);
            
        else
            ihdl2{x} = patch([xvalues fliplr(xvalues)], ribbon1,plot_colors(pc_idx(x),:));
            alpha(0.5);
            hold on;
            % This plots a Line for the Mean
            ihdl1{x} = plot(xvalues,mean_ge,'Color',plot_colors(pc_idx(x),:),'LineWidth',3);
        end
        
    end % for y
    
    % Plot a line going through 0.
    if isempty(strfind(ERPn,'ASSR')) && isempty(strfind(ERPn,'SSVEP'))
        plot([0 0], [ylimits(1) ylimits(2)],'--k');
        plot([xvalues(1) xvalues(end)],[0 0],'--k');
    end
    
    % Add in the tickmarks marking significant time bins.
    % get from sig_times{x}:
    for y = 1:length(sig_times{x})
        plot([xvalues(sig_times{x}(y)) xvalues(sig_times{x}(y))],...
            [ylimits(1)+abs((ylimits(1)*0.1)) ylimits(1)],'LineWidth',3,'Color','k');
    end % for y
    
    %         for y = 1:length(ihdl1)
    %             set(ihdl1{y},'ButtonDownFcn','copyaxis');
    %             set(ihdl2{y},'ButtonDownFcn','copyaxis');
    %         end % for x
    
    set(gca,'XLim',[xvalues(1) xvalues(end)],'YLim',ylimits,'TickDir','out');
    axis([xvalues(1) xvalues(end) ylimits(1) ylimits(2)])
    if length(grand_erp) == 6
        title([ERPn,': Cluster ',num2str(x),' Difference Waves']);
    else
        title([ERPn,': Cluster ',num2str(x)]);
    end
    %legend(binnames);
    
    % Also, plot a blank topo showing which channels were included in the time
    % waveform average.
    cfg = [];
    cfg.layout = ['Biosemi',num2str(num_chans),'.mat'];
    cfg.highlight = 'on';
    cfg.highlightchannel = top_sig_chans{x};
    cfg.highlightcolor = map(color_idx(x),:);
    cfg.style = 'blank';
    cfg.comment = 'no';
    figure
    ft_topoplotER(cfg,tftstruct); % I don't think the 2nd input matters,
    % since it will plot a blank topo, but check.
    colormap jet;
end % for x

if length(grand_erp) == 6 % 6 conditions found
    % make 2 figures, 1 with difference wave topos (done above);
    % and 1 with the original erp data from each condition... do here.
    for x = 1:length(top_sig_chans)
        figure
        ihdl1 = {};
        ihdl2 = {};
        %for y = 1:length(grand_erp)
        for y = 1:4
            % Pull out the relevant channels from grand_erp for this
            % condition and average (for the mean on the graphs)
            temp_ge = grand_erp{y+2}.avg(top_sig_chans{x},:);
            mean_ge = mean(temp_ge,1);
            
            % Also pull out the relevant channels from idata for this
            % condition and put into calculate_erp_variance to get the
            % variance to plot
            temp_idata = idata{y+2}(top_sig_chans{x},:,:);
            temp_iavg = mean(temp_idata,1); % average across the selected channels.
            [var_plot] = calculate_erp_variance({temp_iavg}, 2);
            % flag set to 2... Does across-subjects SEM for variance,
            % calculated on the individual subjects' ERP averaged across
            % the selected channels.
            
            % Now, plot the ribbon for this cluster:
            ribbon1a = mean_ge + var_plot{1};
            ribbon1b = mean_ge - var_plot{1};
            ribbon1 = [ribbon1b fliplr(ribbon1a)];
            
            if length(grand_erp) > 1
                ihdl2{y} = patch([xvalues fliplr(xvalues)], ribbon1,plot_colors(pc_idx(y),:));
                alpha(0.5);
                hold on;
                % This plots a Line for the Mean
                ihdl1{y} = plot(xvalues,mean_ge,'Color',plot_colors(pc_idx(y),:),'LineWidth',3);
                
            else
                ihdl2{x} = patch([xvalues fliplr(xvalues)], ribbon1,plot_colors(pc_idx(x),:));
                alpha(0.5);
                hold on;
                % This plots a Line for the Mean
                ihdl1{x} = plot(xvalues,mean_ge,'Color',plot_colors(pc_idx(x),:),'LineWidth',3);
            end
            
        end % for y
        
        % Plot a line going through 0.
        if isempty(strfind(ERPn,'ASSR')) && isempty(strfind(ERPn,'SSVEP'))
            plot([0 0], [ylimits(1) ylimits(2)],'--k');
            plot([xvalues(1) xvalues(end)],[0 0],'--k');
        end
        
        % Add in the tickmarks marking significant time bins.
        % get from sig_times{x}:
        for y = 1:length(sig_times{x})
            plot([xvalues(sig_times{x}(y)) xvalues(sig_times{x}(y))],...
                [ylimits(1)+abs((ylimits(1)*0.1)) ylimits(1)],'LineWidth',3,'Color','k');
        end % for y
        
        %         for y = 1:length(ihdl1)
        %             set(ihdl1{y},'ButtonDownFcn','copyaxis');
        %             set(ihdl2{y},'ButtonDownFcn','copyaxis');
        %         end % for x
        
        set(gca,'XLim',[xvalues(1) xvalues(end)],'YLim',ylimits,'TickDir','out');
        axis([xvalues(1) xvalues(end) ylimits(1) ylimits(2)])      
        title([ERPn,': Cluster ',num2str(x)]);
        %legend(binnames);
        
%         % Also, plot a blank topo showing which channels were included in the time
%         % waveform average.
%         cfg = [];
%         cfg.layout = ['Biosemi',num2str(num_chans),'.mat'];
%         cfg.highlight = 'on';
%         cfg.highlightchannel = top_sig_chans{x};
%         cfg.highlightcolor = map(color_idx(x),:);
%         cfg.style = 'blank';
%         cfg.comment = 'no';
%         figure
%         ft_topoplotER(cfg,tftstruct); % I don't think the 2nd input matters,
%         % since it will plot a blank topo, but check.
%         colormap jet;
    end % for x            
end % if length(grand_erp) == 6