%% 

function [mean_error, median_error, true_positions, decoded_positions] = decode_by_day1(mouse_folder)

%% Identify files
session_files = get_file_paths_all(mouse_folder); % All files
cell_reg_file_id = contains(session_files, {'cell_reg','cellreg'},'IgnoreCase',1);% Cell reg file index
cell_reg_file = session_files{cell_reg_file_id}; % cell reg file
day1_id = 1;
session_files = session_files(~cell_reg_file_id);
testing_id = 2:numel(session_files);
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
% For each testing session, train a decoder with the session 1 common
% cells, and calculate the decoding errors.
for i =1:numel(testing_id)
    disp(['Calculating decoding errors by Day 1 for mouse ' mouse_id, ' Session ' num2str(testing_id(i)) ' of ' num2str(numel(session_files))])
    [sample_traces_binarized_i, sample_true_position_i, training_traces_binarized_i, training_position_i] = ...
        prepare_bayesian_decode_across_days(mouse_folder, day1_id, testing_id(i), n_spatial_bins, cell_reg_mtx);
    [decoded_position_i, ~, ~] = bayesian_decode_Ca(sample_traces_binarized_i, training_traces_binarized_i, training_position_i);
    mean_error(i) = mean(abs(sample_true_position_i - decoded_position_i));
    median_error(i) = median(abs(sample_true_position_i - decoded_position_i));
    true_positions{i} = sample_true_position_i;
    decoded_positions{i} = decoded_position_i;
end


