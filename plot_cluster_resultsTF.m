function plot_cluster_resultsTF(stat,grandtf,zlimits,p,pt)
% plot_cluster_resultsTF.m
% currently only handles 1 condition...
%
% K. Backer, Feb 2017.

% Check to see if zlimits is a cell array.
% If so, take the largest ranges:
if iscell(zlimits)
    z = cell2mat(zlimits);
    zlimits = [min(z) max(z)];
end

% custom colormap:
% light grey = n.s.
% positive clusters = warm colors
% negative clusters = cool colors
map = [0,1,0.7; 0,1,0; 0,0.5,0; 0,0,1; 0,0.5,1; 0,1,1; 0,0.5,0.5;...
    0.7,0.7,0.7; ...
    0.5,0.5,0; 1,0.5,0; 1,0,0; 0.5,0,0.5; 0.75,0,0.75;1,0,1; 1,1,0];

% Get the indices for the significant positive and negative clusters:
% Similar to ERPs, but now 3D.
color_idx = [];
pos_col_idx = [9:15];
neg_col_idx = fliplr([1:7]);
zero_col_idx = pos_col_idx(1)-1;
rplot = zeros(size(stat.mask));
ppv = 1;
npv = -1;
thold = stat.cfg.alpha;
sc_idxs = {}; % Save the indices for the significant clusters for later plotting.
sc_cnt = 1;
tfbin_thold = 0; % Must have more than 100 t-f bins, where it's significant to be included.
if isfield(stat,'posclusters')
    for x = 1:length(stat.posclusters) % positive clusters.
        if stat.posclusters(x).prob < thold
            % Use stat.posclusterslabelmat to get the cluster indices.
            f = find(stat.posclusterslabelmat==x);
            if length(f)>tfbin_thold
                rplot(f) = ppv;
                if ppv <= length(pos_col_idx)
                    color_idx(sc_cnt) = pos_col_idx(ppv);
                end
                ppv = ppv + 1;
                sc_idxs{sc_cnt} = f;
                sc_cnt = sc_cnt + 1;
            end
        else
            if x == length(stat.posclusters) && ppv == 1;
                display('No positive clusters surpassed p threshold.')
            end
        end % if
    end % for x
else
    display('No positive clusters were found.')
end % if isfield

if isfield(stat,'negclusters')
    for x = 1:length(stat.negclusters) % negative clusters.
        if stat.negclusters(x).prob < thold
            % Use stat.negclusterslabelmat to get the cluster indices.
            f = find(stat.negclusterslabelmat==x);
            if length(f)>tfbin_thold
                rplot(f) = npv;
                if abs(npv) <= length(neg_col_idx)
                    color_idx(sc_cnt) = neg_col_idx(abs(npv));
                end
                npv = npv - 1;
                sc_idxs{sc_cnt} = f;
                sc_cnt = sc_cnt + 1;
            end
        else
            if x == length(stat.negclusters) && npv == -1
                display('No negative clusters surpassed p threshold.')
            end
        end % if
    end % for x
else
    display('No negative clusters were found.')
end % if isfield

num_chans = size(stat.stat,1);
num_freqs = size(stat.stat,2);
num_times = size(stat.stat,3);

% Plot any significant clusters for each channel individually, i.e., in mutliple plots across the scalp: 
if ppv <= length(pos_col_idx) && abs(npv) <= length(neg_col_idx)
    clust_zmap = map;
    clust_zlimits = [-1*((size(map,1)-1)/2) (size(map,1)-1)/2];
    if iscell(pt)
        plot_channel_TFdata({rplot}, num_chans, stat.time, stat.freq,clust_zlimits,[pt{1},' vs ', pt{2}],clust_zmap);
    else
        plot_channel_TFdata({rplot}, num_chans, stat.time, stat.freq,clust_zlimits,pt,clust_zmap);
    end
else
    %clust_zlimits = [-1*((size(map,1)-1)/2) (size(map,1)-1)/2];
    if iscell(pt)
        plot_channel_TFdata({rplot}, num_chans, stat.time, stat.freq,[],[pt{1},' vs ', pt{2}]);
    else
        plot_channel_TFdata({rplot}, num_chans, stat.time, stat.freq,[],pt);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% At this point, use the significant cluster indices to get the time and
% frequencies involved in each cluster.
% because the indices are not in 2D format, make a matrix of channel
% numbers the same size as stat.prob and a matrix of times the same size as
% stat.prob, using repmat.
chan_mat = repmat([1:num_chans]',[1 num_freqs num_times]);
freq_mat = repmat([1:num_freqs],[num_chans 1 num_times]);
% to repmat the 3rd dimension is a bit more complicated...
time_mat = zeros(num_chans, num_freqs,num_times);
for x = 1:num_times
    time_mat(:,:,x) = x;    
end % for x

sig_chans = cell(size(sc_idxs));
sig_times = cell(size(sc_idxs));
sig_freqs = cell(size(sc_idxs));
for x = 1:length(sc_idxs)
    sig_chans{x} = unique(chan_mat(sc_idxs{x}));
    sig_times{x} = unique(time_mat(sc_idxs{x}));
    sig_freqs{x} = unique(freq_mat(sc_idxs{x}));
end % for x
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Another way to try to visualize the cluster results... in one 3D plot...
% Now, manually reshape the stats parameter to a 2D matrix, to
% allow for surface plotting:
% currently, dimensions are: Channel x Frequency x Time
% Want: channel+time x Frequency
% Rows: 1-32 = Ch. 1-32 @ Time Point 1.
% Rows: 33-64 = Ch. 1-32 @ Time Point 2.
% etc.
% Columns = Frequencies 1-N.
rplot2 = zeros(num_chans*num_times, num_freqs);
midx = 1;
% Loop through each time point, extracting the data for all channels and frequencies at that time point.
for t = 1:num_times % Time is 3rd Dimension.
    rplot2(midx:midx+num_chans-1,:) = squeeze(rplot(:,:,t));
    midx = midx+num_chans;
end % for t
% Now, use repmat to make the correct coordinate mappings for each Axis.
Xmat = repmat([1:num_chans]',num_times,num_freqs);
Ymat = zeros(size(Xmat));
yidx = 1;
for y = 1:num_times
    Ymat(yidx:yidx+num_chans-1,1:num_freqs) = repmat(stat.time(y),num_chans,num_freqs);
    yidx = yidx + num_chans;
end
Zmat = repmat(stat.freq,num_chans*num_times,1);
% figure, surf(Xmat,Ymat,Zmat,rplot2) % Whole volume... can't see anything
% but volume faces.

% So, pull out just the significant (non-zero) points, and their corresponding
% coordinates, for plotting.
pidx = find(rplot2);
%figure,scatter3(Xmat(pidx),Ymat(pidx),Zmat(pidx),10,rplot2(pidx));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Now, use stat.mask to essentially threshold grand TF data, to show just
% the activity that is significant.  Do this for each condition, in case
% more than 1 cell of grandtf.
trimmed_data = cell(size(grandtf));
for x = 1:length(grandtf)
    % First, equate size of the TF data and stats (since definitely some time
    % and likely some frequencies were not included in the stats analysis):
    f1 = find(grandtf{x}.freq == stat.freq(1));
    f2 = find(grandtf{x}.freq == stat.freq(end));
    t1 = find(grandtf{x}.time == stat.time(1));
    t2 = find(grandtf{x}.time == stat.time(end));
    trimmed_data{x} = eval(['grandtf{x}.',p,'(:,f1:f2,t1:t2)']);
     
    % Now, mask Trimmed data using stat.mask:
    trimmed_data{x}(find(~stat.mask)) = 0;
    if length(grandtf) == 1 % Plot if only 1 condition.
        if iscell(pt)
            plot_channel_TFdata(trimmed_data(x), num_chans, stat.time, stat.freq, zlimits,pt{x});
        else
            plot_channel_TFdata(trimmed_data(x), num_chans, stat.time, stat.freq, zlimits,pt);
        end
    end % if length(grandtf)
end % for x

% Also, if grandtf contains 2 conditions, maybe it's good to plot the
% difference spectrograms using the trimmed data:
if length(grandtf) == 2
    trimmed_diff = trimmed_data{1} - trimmed_data{2};
    %plot_channel_TFdata({trimmed_diff}, num_chans, stat.time, stat.freq, zlimits, ['Masked: ',pt{1},' - ',pt{2}]);
    
    % Also, plot the raw, unmasked differences between the two:
    d1 = eval(['grandtf{1}.',p,'(:,f1:f2,t1:t2)']);
    d2 = eval(['grandtf{2}.',p,'(:,f1:f2,t1:t2)']);
    %plot_channel_TFdata({d1-d2}, num_chans, stat.time, stat.freq, zlimits, ['Raw: ',pt{1},' - ',pt{2}]);
    
    a = 0.5; % Transparency value.
    plot_channel_TFdataMask([d1-d2], trimmed_diff, num_chans, stat.time, stat.freq, zlimits, [pt{1},' - ',pt{2}],a);
    plot_channel_TFdataMask(d1, trimmed_data{1}, num_chans, stat.time, stat.freq, zlimits, pt{1},a);
    plot_channel_TFdataMask(d2, trimmed_data{2}, num_chans, stat.time, stat.freq, zlimits, pt{2},a);
end % if 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Now, plot the scalp topographies of the TF data for each cluster:
% Set up the general config options for plotting topos of the actual data:
cfg = [];
cfg.parameter = p;
cfg.marker = 'on';
cfg.markersymbol = 'o';
cfg.markersize = 2;
cfg.highlight = 'on';
cfg.highlightsymbol = '*';
cfg.layout = ['Biosemi',num2str(length(stat.label)),'.mat'];
cfg.comment = 'xlim';
cfg.commentpos = 'title';
cfg.colorbar = 'yes';
cfg.zlim = zlimits;

% This one is averaged across all times and all frequencies --> amplitude
% is really washed out...
if length(grandtf) == 1 || length(grandtf) == 2 % Only 1 condition (1-sample t-test)
    % OR 2 conditions found (paired-samples or independent samples t-test)
    sub_cnt = 1;
    figure
    % Loop through each significant cluster:
    for x = 1:length(sc_idxs)
        for y = 1:length(grandtf) % loop through each cell of grandtf (i.e., each condition)
            subplot(length(sc_idxs),length(grandtf),sub_cnt);
            cfg.xlim = [stat.time(sig_times{x}(1)) stat.time(sig_times{x}(end))]; % Times to average across.
            cfg.ylim = [stat.freq(sig_freqs{x}(1)) stat.freq(sig_freqs{x}(end))]; % Freqs to average across.
            cfg.highlightchannel = ft_channelselection(sig_chans{x},grandtf{y}.label);
            ft_topoplotER(cfg,grandtf{y});            
            colormap jet;
            sub_cnt = sub_cnt + 1;
        end % for y
    end % for x    
end % if

% Another way of plotting the cluster topos is to use the rplot matrix as a mask, and
% average across the times and frequencies that are significant, just the
% significant time-frequency bins in each channel.
cidxs = setxor(unique(rplot),0);
if length(grandtf) == 1 || length(grandtf) == 2 % Only 1 condition (1-sample t-test)
    % OR 2 conditions found (paired-samples or independent samples t-test)
    sub_cnt = 1;
    figure
    % Loop through each significant cluster:
    for x = 1:length(cidxs)        
        ch_avgs = cell(size(grandtf));
        for c = 1:length(ch_avgs)
            ch_avgs{c} = zeros(num_chans,1);
        end
        for y = 1:num_chans
            temp_rplot = squeeze(rplot(y,:,:));
            for z = 1:length(grandtf)
                % First, use trimmed data... maybe not all times or all frequencies
                % were inputted into the stats.
                temp_td = squeeze(trimmed_data{z}(y,:,:));
                % Note trimmed data has already been masked, using
                % stat.mask, but it may contain activation from multiple
                % clusters, if more than 1 cluster was significant.
                % So, Mask temp_td with rplot:
                f = find(temp_rplot == cidxs(x));
                ch_avgs{z}(y) = mean(temp_td(f));  
            end % for z
        end % for y
        
        % Now, plot the TOPO:
        for c = 1:length(ch_avgs)
            all_data = grandtf{c};
            if isfield(all_data,'powspctrm')
                all_data = rmfield(all_data,'powspctrm');
            end
            if isfield(all_data,'itpc')
                all_data = rmfield(all_data,'itpc');
            end
            all_data = rmfield(all_data,'freq');
            all_data.time = all_data.time(1);
            all_data.dimord = 'time_chan';
            all_data.avg =ch_avgs{c}';
            subplot(length(cidxs),length(ch_avgs),sub_cnt);
            cfg.xlim = [];
            cfg.ylim = [];
            cfg.comment = 'no';
            cfg.parameter = 'avg';
            ft_topoplotER(cfg,all_data);
            colormap jet;
            sub_cnt = sub_cnt + 1;
        end
    end % for x
    
end % if length

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%