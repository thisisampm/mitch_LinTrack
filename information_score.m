function [IC, PC_idx, place_field_peak, pf_width] = information_score(behavior_mtx, traces, inclusion_vector, direction_preference, varargin)
% Calculate the information score Skaggs et al. 1994 for each neuron

%% Check inputs
if ~isempty(varargin)
    trace_type = varargin{1};
else
    trace_type = 3; % trace type used is C (traces matrix is divided S, C, C_raw)
end

switch trace_type
    case 1
        sampling_rate = 20;
        trace_mtx = traces(:,:,1);
    case 2
        trace_mtx = traces(:,:,2);
    case 3
        sampling_rate = 20;
        z_threshold = 2;
        trace_mtx = extract_binary(traces(:,:,3),sampling_rate, z_threshold);
end

%{
%% Initialize
ncells = size(traces,1); % First dimension of traces is the number of cells
nbins = 40; % Divide track into 40 spatial bins
behavior_mtx = behavior_mtx(inclusion_vector,:);
trace_mtx = trace_mtx(:,inclusion_vector);
%% Extract spatial firing metrics
[activity_by_bin, mean_activity, occupancy] = extract_1D_spatial_tuning(behavior_mtx, trace_mtx, trace_type, nbins);
%}
%traces = traces(:,:,trace_type);

%% Initialize
nbins = 40;
xbinned = discretize(behavior_mtx(inclusion_vector,2),nbins); % Bin the position data
bins = unique(xbinned(~isnan(xbinned))); % Get the bin IDs (a column vector)

% Remove excluded portions
trace_mtx = trace_mtx(:,inclusion_vector);
%xbinned = xbinned(inclusion_vector);

occupancy = sum(xbinned == bins'); % Number of timesteps spent in each bin
px = occupancy/sum(occupancy); % Spatial probability function expressed as fraction of timesteps in each bin

ncells = size(traces,1); % 1st dimension of traces represents number of neurons

lambda_x = nan(ncells,nbins); % Preallocate a matrix to contain average activity of each neuron in each bin


%% Calculate mean Ca activity in each bin for every neuron.
for ibin = bins'
    in_bin = xbinned == ibin; % logical index of current bin occupancy
    if trace_type == 1 % if spikes (S), compute rate
        lambda_x(:,ibin) = sum(trace_mtx(:,in_bin),2)./(sum(in_bin)./sampling_rate); % Rate_in_bin = spikes_in_bin/time_in_bin
    elseif trace_type == 2 % if Ca flouresence, compute mean
        lambda_x(:,ibin) = mean(trace_mtx(:,in_bin),2,'omitnan');
    elseif trace_type == 3 % if C_raw binarized count active frames
        lambda_x(:,ibin) = sum(trace_mtx(:,in_bin),2)./(sum(in_bin));
    else
        error('Trace type not indicated as S (1), C (2), or C_raw (3).');
    end
end

lambda = sum(lambda_x.*px,2); % Calculate the average Ca activity for each neuron across the whole track
log_lambda = log2(lambda_x./lambda); % Calculate log portion of the equation
log_lambda(isinf(log_lambda)) = 0; % log2 of Zero activity bins results in -inf, however, they will be set to zero by 
% also being zero in the first part of equation, so just set to zero now to avoid error resulting from 0*-inf

IC = sum((lambda_x./lambda).*log_lambda.*px,2); % Calculate the information content of each cell's spatial activity

%{
log2_activity = log2(activity_by_bin./mean_activity); % Calculate log portion of the equation
log2_activity(isinf(log2_activity)) = 0; % log2 of Zero activity bins results in -inf, however, they will be set to zero by 
% also being zero in the first part of equation, so just set to zero now to avoid error resulting from 0*-inf

IC = sum((activity_by_bin./mean_activity).*log2_activity.*occupancy,2); % Calculate each cell's spatial information content (IC)
%}
%% Shuffle calcium transient positions and resulting spatial information scores

nshuf = 1000; % Do 1000 shuffles
IC_shuf = nan(ncells, nshuf); % Preallocate neurons x shuffle sized matrix to store IC scores for each cell during each shuffle.


lambda_x_shuf = nan(ncells, nbins, nshuf); % Preallocate a matrix that will store shuffled bin activity during each shuffle

for ishuf = 1:nshuf
    %fprintf('Shuffle %g of %g.\n',i,n_shuf);
    Ca_trial_shuf = circshift(trace_mtx,randi(size(trace_mtx,2)),2); % For each iteration, circularly shuffle the calcium activity to decouple calcium timestamps with position
    for ibin = 1:nbins
        in_bin = xbinned == ibin;
        if trace_type == 1 % if spikes, compute rate
            lambda_x_shuf(:,ibin,ishuf) = sum(Ca_trial_shuf(:,in_bin),2)./(sum(in_bin)/sampling_rate);
        elseif trace_type == 2 % if ca flouresence compute mean
            lambda_x_shuf(:,ibin,ishuf) = mean(Ca_trial_shuf(:,in_bin),2); % For each spatial bin, calculate the average flourescence of each cell 
        elseif trace_type == 3
            lambda_x_shuf(:,ibin,ishuf) = sum(Ca_trial_shuf(:,in_bin),2)./(sum(in_bin)); % If binarized divide active timestamps in bin by total timestamps in bin.
        end
    end
    lambda_shuf = sum(lambda_x_shuf(:,:,ishuf).*px,2);
    log_lambda_shuf = log2(lambda_x_shuf(:,:,ishuf)./lambda_shuf); % Calculate log portion of the information equation
    log_lambda_shuf(isinf(log_lambda_shuf)) = 0; % Set to zero to avoid -inf error, would be mutlipled by a zero anyways.
    IC_shuf(:,ishuf) = sum((lambda_x_shuf(:,:,ishuf)./lambda_shuf).*log_lambda_shuf.*px,2); % Store the information content for each cell for each shuffle iteration
end

%{
shuffled_activity_by_bin = nan(ncells, nbins, nshuf); % Preallocate a matrix to store resulting tuning curve on each shuffle

for ishuf = 1:nshuf
    trace_mtx_shuffled = circshift(trace_mtx,randi(size(trace_mtx,2)),2); % For each iteration, circularly shuffle the calcium activity to decouple calcium timestamps with position
    [shuffled_activity_by_bin(:,:,ishuf), ~, ~] = extract_1D_spatial_tuning(behavior_mtx, trace_mtx_shuffled, trace_type, nbins); % Calculate shuffled tuning curves
    log2_shuffled_activity = log2(shuffled_activity_by_bin(:,:,ishuf)./mean_activity); % Calculate log portion of the equation
    log2_shuffled_activity(isinf(log2_shuffled_activity)) = 0;
    IC_shuf(:,ishuf) = sum((shuffled_activity_by_bin(:,:,ishuf).*mean_activity).*log2_shuffled_activity.*occupancy,2);
end
%}
sig_SI_idx = IC > prctile(IC_shuf,95,2); % An index of cells with greater spatial information than 95% of information scores resulting from shuffles



%% Determine if a cell is active on a minimum number of trials
[trials, ~, direction] = get_trials(behavior_mtx); % All trials for next function

switch direction_preference
    case 'none'
        trials = trials(direction == 0 | direction == 1);
    case 'left'
        trials = trials(direction == 1);
    case 'right'
        trials = trials(direction == 0);
end
        
        
rate_mtx_3d = binned_trl_activity_MDS(behavior_mtx,traces,trace_type, nbins, trials); % Generate a trial x bin x neuron matrix

min_active_trials = round(0.3*numel(trials)); % set the minimum number of trials for place cell consideration
cell_active_idx = false(size(traces,1),1); % Preallocate an index for active threshold criteria

for i_cell = 1:size(rate_mtx_3d,3)
    cell_active = rate_mtx_3d(:,:,i_cell) ~= 0; % Find active bins
    n_trial_active = max(sum(cell_active)); % Find how many trials there are activity
    if n_trial_active > min_active_trials
        cell_active_idx(i_cell) = 1; % Indicate if cell was active on more than threshold trials
    end
end

%% Determine the width of the place field
pvalues = sum(lambda_x_shuf > lambda_x,3)/nshuf; % p value determines where tuning curve is significantly above shuffled tuning distribution
%pvalues = sum(shuffled_activity_by_bin > activity_by_bin,3)/nshuf; % p value determines where tuning curve is significantly above shuffled tuning distribution
sig_pval = pvalues < 0.05; % Determine where tuning curves are over 
sig_pval = [zeros(ncells,1) sig_pval zeros(ncells,1)]; % Bookend with zeros for use with diff

% For each cell with activity define place field with as the maximum consecutive above 
pf_width = nan(ncells,1);
active_cells = find(sig_SI_idx & cell_active_idx); % Place cells
for i = active_cells'
    pf_start = find(diff(sig_pval(i,:)) == 1); % Determine bins where the pvalue goes from above threshold to below it
    pf_end = find(diff(sig_pval(i,:)) == -1); % Determine bins where pvalue goes from under threshold to above
    pf_width(i) = max(pf_end - pf_start); % Take the place field width as the largest consecutive section below p threshold
end

        
%% Output variables
PC_idx = (sig_SI_idx & cell_active_idx); % Report place cells as those with sig SI and activity on a minimum number of trials
[~ ,place_field_peak] = max(lambda_x,[],2); % Report the placement of place field centroid





