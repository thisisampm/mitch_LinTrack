function [IC, PC_idx, place_field_peak] = information_score2(behavior_mtx, traces, varargin)
% Calculate the information score Skaggs et al. 1994 for each neuron

%% Check inputs
if ~isempty(varargin)
    trace_type = varargin{1};
else
    trace_type = 2; % trace type used is C (traces matrix is divided S, C, C_raw)
end

%traces = traces(:,:,trace_type);

if trace_type == 1
    sampling_rate = 100; % If using spikes (S), define the sampling rate. Here it is 100 Hz because 
    % that's what it was resampled to when aligning and interpolating
    % traces with behaviour data.
    fprintf('Information score metrics calculated using spikes (S).\n')
elseif trace_type == 2
    fprintf('Information score metrics calculated using C df/F.\n')
elseif trace_type == 3
    fprintf('Information score metrics calculated using binarized C_raw.\n')
    sampling_frequency = 20;
    z_threshold = 2;
    traces(:,:,3) = extract_binary(traces(:,:,trace_type),sampling_frequency,z_threshold);
end

%% Initialize
%{ 
%Update analysis to analyze each direction 
in_trial_idx = ~isnan(behavior_mtx(:,4)); % determine when the mouse is running a trial
nbins = 40; % 40 spatial bins across the track
xbinned = discretize(behavior_mtx(in_trial_idx,2),nbins); % Divide position during trials into 40 bins
%xpos_binned = discretize(behavior_mtx(:,2),nbins); % Divide position on track into 40 bins
remove_xpos_nans = ~isnan(xbinned);
xbinned = xbinned(remove_xpos_nans);
px = occupancy/numel(xbinned); % Calculate spatial probability distribution as the fraction of time spent in each bin
ca_during_trial = traces(:,in_trial_idx,trace_type);
lambda_x = nan(size(ca_during_trial, 1), nbins); % Preallocate a neurons x bins matrix to store each cell's mean Ca activity in each spatial bin
%}

speed_threshold = 5; % Set a 5 cm/s cutoff
running_ts = behavior_mtx(:,5) > speed_threshold; % Determine when the mouse is running above threshold
leftward = isolate_direction(behavior_mtx(:,2),'left'); % Determine leftward direction
rightward = isolate_direction(behavior_mtx(:,2),'right'); % Repeat for rightward
xbinned = discretize(behavior_mtx(:,2),40); % Bin the position data
bins = unique(xbinned(~isnan(xbinned))); % Get the bin IDs (a column vector)right_occpancy = sum

left_occupancy = sum(xbinned(running_ts & leftward) == bins'); % Determine frames spent in each spatial bin during trial
right_occupancy = sum(xbinned(running_ts & rightward) == bins'); % Determine frames spent in each spatial bin during trial


%Ca_during_trial = traces(:,remove_xpos_nans,3);

%% Calculate mean Ca activity in each bin for every neuron.
for ibin = unique(xbinned)'
    in_bin = xbinned == ibin;
    if trace_type == 1 % if spikes (S), compute rate
        lambda_x(:,ibin) = sum(ca_during_trial(:,in_bin),2)./(sum(in_bin)/sampling_rate); % Rate_in_bin = spikes_in_bin/time_in_bin
        % time_in_bin = frames_in_bin/frame_rate
    elseif trace_type == 2 % if Ca flouresence, compute mean
        lambda_x(:,ibin) = mean(ca_during_trial(:,in_bin),2);
    elseif trace_type == 3 % if C_raw binarize and count
        lambda_x(:,ibin) = sum(ca_during_trial(:,in_bin),2)./(sum(in_bin));
    else
        error('Trace type not indicated as S (1), C (2), or C_raw (3).');
    end
end

lambda = sum(lambda_x.*px,2); % Calculate the average Ca activity for each neuron across the whole track
log_lambda = log2(lambda_x./lambda); % Calculate log portion of the equation
log_lambda(isinf(log_lambda)) = 0; % log2 of Zero activity bins results in -inf, however, they will be set to zero by 
% also being zero in the first part of equation, so just set to zero now to avoid error resulting from 0*-inf

IC = sum((lambda_x./lambda).*log_lambda.*px,2); % Calculate the information content of each cell's spatial activity

%% Shuffle calcium transient positions and resulting spatial information scores

n_shuf = 1000; % Do 500 shuffles
IC_shuf = nan(size(ca_during_trial,1),n_shuf); % Preallocate neurons x shuffle sized matrix to store IC scores for each cell during each shuffle.

lambda_x_shuf = nan(size(ca_during_trial,1),nbins); % Preallocate a matrix that will store shuffled bin activity during each shuffle

for i = 1:n_shuf
    %fprintf('Shuffle %g of %g.\n',i,n_shuf);
    Ca_trial_shuf = circshift(ca_during_trial,randi(size(ca_during_trial,2)),2); % For each iteration, circularly shuffle the calcium activity to decouple calcium timestamps with position
    for ibin = unique(xbinned)'
        in_bin = xbinned == ibin;
        if trace_type == 1 % if spikes, compute rate
            lambda_x_shuf(:,ibin) = sum(Ca_trial_shuf(:,in_bin),2)./(sum(in_bin)/sampling_rate);
        elseif trace_type == 2 % if ca flouresence compute mean
            lambda_x_shuf(:,ibin) = mean(Ca_trial_shuf(:,in_bin),2); % For each spatial bin, calculate the average flourescence of each cell 
        elseif trace_type == 3
            lambda_x_shuf(:,ibin) = sum(Ca_trial_shuf(:,in_bin),2)./(sum(in_bin)); % If binarized divide active timestamps in bin by total timestamps in bin.
        end
    end
    lambda_shuf = sum(lambda_x_shuf.*px,2);
    log_lambda_shuf = log2(lambda_x_shuf./lambda_shuf); % Calculate log portion of the information equation
    log_lambda_shuf(isinf(log_lambda_shuf)) = 0; % Set to zero to avoid -inf error, would be mutlipled by a zero anyways.
    IC_shuf(:,i) = sum((lambda_x_shuf./lambda_shuf).*log_lambda_shuf.*px,2); % Store the information content for each cell for each shuffle iteration
end

sig_SI_idx = IC > prctile(IC_shuf,95,2); % An index of cells with greater spatial information than 95% of information scores resulting from shuffles

%% Determine if a cell is active on a minimum number of trials
trials = unique(behavior_mtx(~isnan(behavior_mtx(:,4)),4)); % All trials for next function
rate_mtx_3d = binned_trl_activity_MDS(behavior_mtx,traces,trace_type, nbins, trials); % Generate a trial x bin x neuron matrix
min_active_trials = 3; % set the minimum number of trials for place cell consideration
cell_active_idx = false(size(traces,1),1); % Preallocate an index for active threshold criteria

for i_cell = 1:size(rate_mtx_3d,3)
    cell_active = rate_mtx_3d(:,:,i_cell) ~= 0; % Find active bins
    n_trial_active = max(sum(cell_active)); % Find how many trials there are activity
    if n_trial_active > min_active_trials
        cell_active_idx(i_cell) = 1; % Indicate if cell was active on more than threshold trials
    end
end
        
%% Output variables
PC_idx = (sig_SI_idx & cell_active_idx); % Report place cells as those with sig SI and activity on a minimum number of trials
[~ ,place_field_peak] = max(lambda_x,[],2); % Report the placement of place field centroid





