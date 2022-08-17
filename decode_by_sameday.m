function [mean_error, median_error, true_positions, decoded_positions] = decode_by_sameday(mouse_folder)

%% Identify files
session_files = get_file_paths_all(mouse_folder); % All files
cell_reg_file_id = contains(session_files, {'cell_reg','cellreg'},'IgnoreCase',1);% Cell reg file index
cell_reg_file = session_files{cell_reg_file_id}; % cell reg file

session_files = session_files(~cell_reg_file_id);
testing_id = 1:numel(session_files);
[~, mouse_id] = fileparts(mouse_folder); % Get the mouse string name for later

load(cell_reg_file);

if ~exist('cell_regist_mtx','var')
    cell_reg_mtx = cell_registered_struct.cell_to_index_map;
elseif exist('cell_regist_mtx','var')
    cell_reg_mtx = cell_regist_mtx;
end

%% Train decoder and get errors for each day
% Preallocate for storage
mean_error = nan(numel(testing_id),1);
median_error = nan(numel(testing_id),1);
true_positions = cell(numel(testing_id),1);
decoded_positions = cell(numel(testing_id),1);

n_spatial_bins = 40;
training_fraction = 0.7; % Fraction of trials to use for training

% For each testing session, train a decoder with a fraction of the trials
% and test on the remaining fraction
for i =1:numel(testing_id)
    disp(['Calculating same day decoding error ' mouse_id, ' Session ' num2str(testing_id(i)) ' of ' num2str(numel(session_files))])
    load(session_files{i},'behavior_mtx','traces');
    session_nans = (logical(sum(isnan(traces(:,:,3)))) | isnan(behavior_mtx(:,2))'); % identify NaNs that will block binarization
    behavior_mtx = behavior_mtx(~session_nans,:);
    traces = traces(:,~session_nans,:);
    [sample_traces_binarized_i, sample_true_position_i, training_traces_binarized_i, training_position_i] = ...
        prepare_bayesian_decode(behavior_mtx, traces(:,:,3), training_fraction, n_spatial_bins);
    [decoded_position_i, ~, ~] = bayesian_decode_Ca(sample_traces_binarized_i, training_traces_binarized_i, training_position_i);
    mean_error(i) = mean(abs(sample_true_position_i - decoded_position_i));
    median_error(i) = median(abs(sample_true_position_i - decoded_position_i));
    true_positions{i} = sample_true_position_i;
    decoded_positions{i} = decoded_position_i;
end
