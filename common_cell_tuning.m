function tuning_curve_matrix = common_cell_tuning(session_files, cell_regist_mtx, sort_by_sesh)
% plot the tuning curves of cells that are common across all session files
% sort rows by peak with reference to session number sort_by_sesh

%for session files, try: get_file_paths_targeted('C:\Users\ampm1\Documents\MATLAB\mitch_LinTrack\control\152-2', {'.mat'})

% details
num_spatial_bins = 40;
traces_type = 2; % C

% input
num_sessions = size(session_files,1);

% common cells
common_cells = sum(cell_regist_mtx>0,2)==num_sessions;

% iterate through sessions computing tuning curves for all cells
tuning_curve_matrix = nan(sum(common_cells), num_spatial_bins, num_sessions);
for isesh = 1:num_sessions
    
    % print update
    disp(['Starting session ' num2str(isesh)])
    
    % load session
    load(session_files{isesh}, 'behavior_mtx', 'traces')
    
    % load tuning curves
    [~, unsorted_session_tuning_curves] = sort_cell_activity(behavior_mtx, traces(:,:,traces_type), num_spatial_bins);
    
    % common cells only
    common_cell_ids = cell_regist_mtx(common_cells,isesh);
    tuning_curve_matrix(:,:,isesh) = unsorted_session_tuning_curves(common_cell_ids,:);
   
    % update figures
    drawnow
end

% sort by peak
[~,sort_idx] = sort_rows_by_peak(tuning_curve_matrix(:,:,sort_by_sesh));
tuning_curve_matrix = tuning_curve_matrix(sort_idx,:,:);

% plot
figure
for isesh = 1:num_sessions
   subplot(1,num_sessions,isesh)
   imagesc(zscore_mtx(tuning_curve_matrix(:,:,isesh)')')
   caxis([-1.5 1.5])
end