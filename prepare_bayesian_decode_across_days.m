

function [sample_traces_binarized, sample_true_position, training_traces_binarized, training_position] ...
    = prepare_bayesian_decode_across_days(mouse_folder, training_session_id, testing_session_id, n_spatial_bins, varargin)

%%%%%%%%% INPUTS %%%%%%%%%%%
% mouse_folder: a string leading to a directory containing session
% behaviour and calcium data
% training_session_id: which session to use as the training day.
% Sessions to compare: a vector of sessions with the indicating
% which to test
% N spatial bins: the number of spatial bins to divide the track into
% varagin: Optional cell_reg_mtx

%%%%%%%%% OUTPUTS %%%%%%%%%%%
% Sample_traces_binarized: the testing data binarized as active or not
% based on Etter et al. 2020
% Sample_true_position: the actual position for error calculation
% training_traces_binarized: the binarized training data
% training_position: corresponding track position for training

%% Load data
session_files = get_file_paths_all(mouse_folder); % Get the files from the mouse folder
cell_reg_file_id = contains(session_files,{'cellreg','cell_reg'},'IgnoreCase',1); % Index the cell reg file
cell_reg_file = session_files{cell_reg_file_id};% Cell Reg filepath
session_files = session_files(~cell_reg_file_id);
% Identify the session for training and the session for testing
training_file = session_files{training_session_id}; % 
testing_file = session_files{testing_session_id};
% Load the stored data into structures with the fields "behavior_mtx and
% traces"
training_data = load(training_file,'behavior_mtx','traces');
testing_data = load(testing_file,'behavior_mtx','traces');
% Check for NaNs that will be problematic for filtfilt in extract_binary
training_data_nans = (logical(sum(isnan(training_data.traces(:,:,3))))|isnan(training_data.behavior_mtx(:,2))');
training_data.behavior_mtx = training_data.behavior_mtx(~training_data_nans,:);
training_data.traces = training_data.traces(:,~training_data_nans,:);

testing_data_nans = (logical(sum(isnan(testing_data.traces(:,:,3))))|isnan(testing_data.behavior_mtx(:,2))');
testing_data.behavior_mtx = testing_data.behavior_mtx(~testing_data_nans,:);
testing_data.traces = testing_data.traces(:,~testing_data_nans,:);

% If cell reg matrix not given, load the cell reg file for identifying cells across sessions
if ~isempty(varargin)
    cell_reg_mtx = varargin{1};
else
    load(cell_reg_file,"cell_registered_struct")
    cell_reg_mtx = cell_registered_struct.cell_to_index_map;
end

%% Identify trials in both sessions and extract trial position data
sampling_frequency = 20; % Resampled to 100 Hz when loading by AMPM code, 20 Hz with reload

[training_trials, training_trial_markers, ~] = get_trials(training_data.behavior_mtx);
[testing_trials, testing_trial_markers, ~] = get_trials(training_data.behavior_mtx);

% Restrict to trials with sufficient running speed
training_trial_durations = (training_trial_markers(:,2)-training_trial_markers(:,1))/sampling_frequency; % Resampled to 100 Hz, so divide number of frames by 100 to get time in seconds
training_too_long = training_trial_durations > 3; % Index to exclude long trials (longer than 3 s)
training_trials = training_trials(~training_too_long); % Exclude long trials

% Repeat restrictions for testing
testing_trial_durations = (testing_trial_markers(:,2)-testing_trial_markers(:,1))/sampling_frequency; % Resampled to 100 Hz, so divide number of frames by 100 to get time in seconds
testing_too_long = testing_trial_durations > 3; % Index to exclude long trials (longer than 3 s)
testing_trials = testing_trials(~testing_too_long); % Exclude long trials

% Logical index for when running included trials
in_training_trial = logical(sum((training_data.behavior_mtx(:,4) == training_trials'),2)); % Index designating training trial timestamps
in_testing_trial = logical(sum((testing_data.behavior_mtx(:,4) == testing_trials'),2)); % Index testing trial timestamps

% Extract the training and testing position data
training_position = discretize(training_data.behavior_mtx(in_training_trial,2),n_spatial_bins);
sample_true_position = discretize(testing_data.behavior_mtx(in_testing_trial,2), n_spatial_bins);

%% Determine the cells common to both sessions and extract their Ca data
common_cell_index = (cell_reg_mtx(:,training_session_id) > 0 & cell_reg_mtx(:,testing_session_id) > 0); % Logical index for cells active in both sessions
training_cells = cell_reg_mtx(common_cell_index, training_session_id);
testing_cells = cell_reg_mtx(common_cell_index, testing_session_id);

% Binarize the Ca data using Etter et al. 2020 protocol
z_threshold = 2; % Threshold of increase in Ca flour to be registered.
training_traces_binarized_all = extract_binary(zscore_mtx(training_data.traces(training_cells,:,3)')',sampling_frequency,z_threshold);
testing_traces_binarized_all  = extract_binary(zscore_mtx(testing_data.traces(testing_cells,:,3)')',sampling_frequency,z_threshold);

% Indexing the binarized trace data down only to trial running times
training_traces_binarized = training_traces_binarized_all(:,in_training_trial)';
sample_traces_binarized = testing_traces_binarized_all(:,in_testing_trial)';



