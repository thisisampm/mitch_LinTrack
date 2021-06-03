
%% behavior code

% first direction of travel on every trial in a session
[dir_vect] = trial_direction(behavior_mtx);



%% single cell plots

% dot plot
dot_plot(behavior_mtx, traces, neuron);

% activity from each trial binned
[trl_act_mtx] = binned_trl_activity(behavior_mtx, trace_col_vect, num_bins);

% pairwise trial activity (for finding place cells)
[mean_corr, corr_matrix] = pairwise_trial_corrs(behavior_mtx, trial_rate_mtx);
[max_mean_corrs, mean_corrs] = ALL_pairwise_trial_corrs(behavior_mtx, rate_mtx_3d);


%% session-level plots

% all cells binned trial activity matrices
[rate_mtx_3d] = sesh_bin_trl_act(behavior_mtx, trace_mtx, num_bins);

% sort activity by peak
[sorted_rate_mtx, rate_mtx] = sort_cell_activity(behavior_mtx, trace_mtx, num_bins);


%% multi session-level plots

% ziv style place cell plots over sessions
tuning_curve_matrix = common_cell_tuning(session_files, cell_reg_mtx, reference_session);
tuning_curve_matrix = all_cell_tuning(session_files, cell_reg_mtx, reference_session);

% bin to bin spatial correlations for each cell, plotted to show session to
% session averages (z scored raw traces seems to work best)
sesh_to_sesh_corrs(tuning_curve_matrix, reference_session)
% [cor_mtx, ebp_ref, ebp_adj] = sesh_to_sesh_corrs(tc_mtx_3d_raw, 1);


%% cell reg plots

% plot common cells between all sessions
common_cells %script

% plot footprints
trace_footprint
trace_footprints_sessions
trace_footprints_sessions_overlap(spatial_footprints_corrected, cell_regist_mtx)

%% loading data