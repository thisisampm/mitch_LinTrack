function [activity_by_bin, mean_activity, px] = extract_1D_spatial_tuning(behavior_mtx, trace_mtx, trace_type, nbins)


%% Initialize

switch trace_type
    case 1
        sampling_rate = 20;
end

xbinned = discretize(behavior_mtx(:,2),nbins); % Bin the position data
bins = unique(xbinned(~isnan(xbinned))); % Get the bin IDs (a column vector)


occupancy = sum(xbinned == bins'); % Number of timesteps spent in each bin
px = occupancy/sum(occupancy); % Spatial probability function expressed as fraction of timesteps in each bin

ncells = size(trace_mtx,1); % 1st dimension of traces represents number of neurons

activity_by_bin = nan(ncells,nbins); % Preallocate a matrix to contain a likelihood of activity for each neuron given occupancy in each bin

%% Calculate mean Ca activity in each bin for every neuron.
for ibin = bins'
    in_bin = xbinned == ibin; % logical index of current bin occupancy
    if trace_type == 1 % if spikes (S), compute rate
        activity_by_bin(:,ibin) = sum(trace_mtx(:,in_bin),2)./(sum(in_bin)./sampling_rate); % Rate_in_bin = spikes_in_bin/time_in_bin
    elseif trace_type == 2 % if Ca flouresence, compute mean
        activity_by_bin(:,ibin) = mean(trace_mtx(:,in_bin),2,'omitnan');
    elseif trace_type == 3 % if C_raw binarized count active frames
        activity_by_bin(:,ibin) = sum(trace_mtx(:,in_bin),2)./(sum(in_bin));
    else
        error('Trace type not indicated as S (1), C (2), or C_raw (3).');
    end
end

mean_activity = sum(activity_by_bin.*px,2); % Calculate the average Ca activity for each neuron across the whole track


