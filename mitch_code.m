
%% behavior code

% first direction of travel on every trial in a session
[dir_vect] = trial_direction(behavior_mtx);



%% single cell plots

% dot plot
dot_plot(behavior_mtx, traces, neuron);

% activity from each trial binned
[trl_act_mtx] = binned_trl_activity(behavior_mtx, trace_col_vect, num_bins);

% pairwise trial activity (for finding place cells)
trial_rate_mtx = sesh_bin_trl_act(behavior_mtx, traces(:,:,2), 40, unique(behavior_mtx(~isnan(behavior_mtx(:,4)),4)));
[max_mean_corrs, mean_corrs] = ALL_pairwise_trial_corrs(behavior_mtx, trial_rate_mtx); %THIS ONE
[mean_corr, corr_matrix] = pairwise_trial_corrs(behavior_mtx, trial_rate_mtx);

%% session-level plots

% all cells binned trial activity matrices
[rate_mtx_3d] = sesh_bin_trl_act(behavior_mtx, trace_mtx, num_bins);

% sort activity by peak
[sorted_rate_mtx, rate_mtx] = sort_cell_activity(behavior_mtx, trace_mtx, num_bins);


%% multi session-level plots

% ziv style place cell plots over sessions
tuning_curve_matrix = common_cell_tuning(session_files, cell_reg_mtx, reference_session);
tuning_curve_matrix = all_cell_tuning(session_files, cell_reg_mtx, reference_session);
% all_tcm = all_cell_tuning_multi([ {'152-2'} {'159-2'} {'167-4'}], 'control', 1);


% load('cell_regist_152-2.mat', 'cell_regist_mtx')
% sesh_nums = 1:8; 
% all_sesh = get_file_paths_targeted('C:\Users\ampm1\Documents\MATLAB\mitch_LinTrack\data\control\152-2', {'.mat'}); 
%  % or load('tuning_curve_matrix_152-2.mat', 'tc_mtx_3d')
% tc_mtx_3d = all_cell_tuning(all_sesh(sesh_nums), cell_regist_mtx(:,sesh_nums), 1);

% bin to bin spatial correlations for each cell, plotted to show session to
% session averages (z scored raw traces seem to work best)
sesh_to_sesh_corrs(tuning_curve_matrix, reference_session)
% [cor_mtx, ebp_ref, ebp_adj] = sesh_to_sesh_corrs(tc_mtx_3d_raw, 1);


%% multi subject plots
all_tcm = all_cell_tuning_multi([ {'152-2'} {'159-2'} ], 'control', 1);
sesh_to_sesh_corrs_multi(subj_cell, reference_session);

%% cell reg plots

% plot common cells between all sessions
common_cells %script

% plot footprints
trace_footprint
trace_footprints_sessions
trace_footprints_sessions_overlap(spatial_footprints_corrected, cell_regist_mtx)


%% loading data

% finds data located in folder_path and saves it in the mitch_linTrack 
% matlab folder according to mouse and day
load_mitchfiles(folderpath);

[behavior_mtx, traces] = load_mitchdata(foldername); %loads data from a mitch-supplied CNMFE and deepLab cut output folder

[traces, frame_times] = load_cnmfe(fp); % load flor traces and timestamps from cnmfe output