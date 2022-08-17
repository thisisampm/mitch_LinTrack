%% GCamP7f animal driver script

%% Get session files

G712_files = get_file_paths_targeted('C:\Users\mitch\OneDrive - University of Toronto\MATLAB\mitch_LinTrack\irradiation\G7.1-2','2021-');
load('C:\Users\mitch\OneDrive - University of Toronto\MATLAB\mitch_LinTrack\irradiation\G7.1-2\cellRegistered_20211117_170437.mat');
G712_ind = cell_registered_struct.cell_to_index_map;

G714_files = get_file_paths_targeted('C:\Users\mitch\OneDrive - University of Toronto\MATLAB\mitch_LinTrack\control\G7.1-4', '2021-');
G714_files = G714_files(1:end-1);
load('C:\Users\mitch\OneDrive - University of Toronto\MATLAB\mitch_LinTrack\control\G7.1-4\cellRegistered_20211113_095102_0930-1014sessions.mat');
G714_ind = cell_registered_struct.cell_to_index_map;

G721_files = get_file_paths_targeted('C:\Users\mitch\OneDrive - University of Toronto\MATLAB\mitch_LinTrack\exercise\G7.2-1', '2021-');
load('C:\Users\mitch\OneDrive - University of Toronto\MATLAB\mitch_LinTrack\exercise\G7.2-1\cellRegistered_20211122_071728_0930-1014.mat');
G721_ind = cell_registered_struct.cell_to_index_map;

G724_files = get_file_paths_targeted('C:\Users\mitch\OneDrive - University of Toronto\MATLAB\mitch_LinTrack\exercise\G7.2-4', '2021-');
G724_files = G724_files(1:end-1); % drop cell reg file here
load('C:\Users\mitch\OneDrive - University of Toronto\MATLAB\mitch_LinTrack\exercise\G7.2-4\cellRegistered_2021-10-12_low_cell_N_session_dropped.mat')
G724_ind = cell_registered_struct.cell_to_index_map;

%% Run irradiation animal G7.1-2
[tuning_mtx_ref_G712, tuning_mtx_all_G712] = all_cell_tuning(G712_files, G712_ind,1);
[cor_mtx_G712, ebp_reference_G712, ebp_adjacent_G712] = sesh_to_sesh_corrs(tuning_mtx_all_G712, 1);

%% Run control animal G7.1-4
[tuning_mtx_ref_G714, tuning_mtx_all_G714] = all_cell_tuning(G714_files, G714_ind,1);
[cor_mtx_G714, ebp_reference_G714, ebp_adjacent_G714] = sesh_to_sesh_corrs(tuning_mtx_all_G714, 1);

%% Run exercise animal G7.2-1
[tuning_mtx_ref_G721, tuning_mtx_all_G721] = all_cell_tuning(G721_files, G721_ind,1);
[cor_mtx_G721, ebp_reference_G721, ebp_adjacent_G721] = sesh_to_sesh_corrs(tuning_mtx_all_G721, 1);

%% Run exercise animal G7.2-4
[tuning_mtx_ref_G724, tuning_mtx_all_G724] = all_cell_tuning(G724_files, G724_ind,1);
[cor_mtx_G724, ebp_reference_G724, ebp_adjacent_G724] = sesh_to_sesh_corrs(tuning_mtx_all_G724, 1);