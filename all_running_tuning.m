


function [tuning_curves, px] = all_running_tuning(behavior_mtx,trace_mtx, inclusion_vector)
% Use inclusion vector to specify times when the mouse is running and/or in
% a certain direction
% Trace_mtx is a neuron x timestep matrix of binarized activity

nbins = 40;
xbinned = discretize(behavior_mtx(:,2),nbins); % Bin the position data
bins = unique(xbinned(~isnan(xbinned))); % Get the bin IDs (a column vector)

% Remove excluded portions
trace_mtx = trace_mtx(:,inclusion_vector);
xbinned = xbinned(inclusion_vector);

occupancy = sum(xbinned == bins'); % Number of timesteps spent in each bin
px = occupancy/sum(occupancy); % Spatial probability function expressed as fraction of timesteps in each bin

ncells = size(trace_mtx,1); % 1st dimension of traces represents number of neurons

tuning_curves = nan(ncells,nbins); % Preallocate a matrix to contain average activity of each neuron in each bin

for ibin = bins'
    in_bin = xbinned == ibin; % logical index of current bin occupancy
    tuning_curves(:,ibin) = sum(trace_mtx(:,in_bin),2)./(sum(in_bin));
end

