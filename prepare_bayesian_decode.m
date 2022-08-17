% prepare_bayesian_decode.m Function

% Function to reshape behaviour and trace data into sample (testing) and
% training datasets for Bayesian decoding based on the proportion of trials
% choosen to be used in testing

% Inputs: 
% behavior_mtx = full 6 column behaviour mtx 
% traces = 2D neuron x timestep matrix with the Ca trace data from entire
% session
% training_fraction = the proportion of trials to be used for training
% n_sample_bins = the number of discret bins to divide the spatial position
% into

function [sample_traces_binarized, sample_true_position, training_traces_binarized, training_position] = prepare_bayesian_decode(behavior_mtx, traces_raw, training_fraction, n_spatial_bins)

%% Identify trials
sampling_frequency = 20; % Resampled to 100 Hz when loading by AMPM, 20 Hz during reload

[trials, trial_markers, direction] = get_trials(behavior_mtx);

% Restrict to trials with sufficient running speed
trial_durations = (trial_markers(:,2)-trial_markers(:,1))/sampling_frequency; % Resampled to 100 Hz, so divide number of frames by 100 to get time in seconds
too_long = trial_durations > 5; % Index to exclude long trials (longer than 3 s)
trials = trials(~too_long); % Exclude long trials
trials = trials(direction(~too_long) == 1); % Restrict direction rightward trials

ntrials = numel(trials); % The number of trials included in a session
training_trials = trials(sort(randperm(ntrials, floor(ntrials*training_fraction)))); % Randomly select training trials based on the predefined fraction
testing_trials = trials(~ismember(trials,training_trials)); % Designate the remaining trials for testing

% Quickly index all trials by creating a matrix of logical vectors and
% collapsing them by summing
in_training_trial = logical(sum((behavior_mtx(:,4) == training_trials'),2)); % Index designating training trial timestamps
in_testing_trial = logical(sum((behavior_mtx(:,4) == testing_trials'),2)); % Index testing trial timestamps

%% Index out relevant trace and position data from the trials
% Use Etter et al. 2020 binarization method. 
z_threshold = 2; % Threshold of increase in Ca flour to be registered.
binarized_traces = extract_binary(zscore_mtx(traces_raw')', sampling_frequency, z_threshold);

binned_position = discretize(behavior_mtx(:,2),n_spatial_bins);

sample_traces_binarized = binarized_traces(:,in_testing_trial)'; % Trace data during sample (testing) trials, output with cells as columns
% sample_true_position = discretize(behavior_mtx(in_testing_trial,2),n_spatial_bins); % True positions to check decoding errors
sample_true_position = binned_position(in_testing_trial);

training_traces_binarized = binarized_traces(:, in_training_trial)'; % Ca trace data during training trials
% training_position = discretize(behavior_mtx(in_training_trial,2),n_spatial_bins); % Corresponding spatial positions during training trials
training_position = binned_position(in_training_trial);

