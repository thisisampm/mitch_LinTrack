function [tuning_curve_matrix_ref, tuning_curve_mtx_all] = all_cell_tuning(session_files, cell_regist_mtx, reference_session)
% plot the tuning curves of cells that are common across all session files
% sort rows by peak with reference to session number sort_by_sesh

%for session files, try: get_file_paths_targeted('C:\Users\ampm1\Documents\MATLAB\mitch_LinTrack\control\152-2', {'.mat'})

% details
num_spatial_bins = 40;
traces_type = 3; % S, C, RAW

% input
num_sessions = size(session_files,1);

% iterate through sessions computing tuning curves for all cells
tuning_curve_mtx_all = nan(size(cell_regist_mtx,1), num_spatial_bins, num_sessions);
for isesh = 1:num_sessions
    
    % print update
    disp(['Starting session ' num2str(isesh)])
    
    % load session
    load(session_files{isesh}, 'behavior_mtx', 'traces')
    
    % load tuning curves
    [~, unsorted_session_tuning_curves] = sort_cell_activity(behavior_mtx, traces(:,:,traces_type), num_spatial_bins);
    
    % common cells only
    tuning_curve_mtx_all(cell_regist_mtx(:,isesh)>0,:,isesh) = unsorted_session_tuning_curves(cell_regist_mtx(cell_regist_mtx(:,isesh)>0,isesh),:);
   
    % update figures
    drawnow
end

% only include cells active during reference session
tuning_curve_matrix_ref = tuning_curve_mtx_all(~isnan(tuning_curve_mtx_all(:,1,reference_session)),:,:);

% sort by peak
[~,sort_idx] = sort_rows_by_peak(tuning_curve_matrix_ref(:,:,reference_session));
tuning_curve_matrix_ref = tuning_curve_matrix_ref(sort_idx,:,:);

% plot
figure
for isesh = 1:num_sessions
   subplot(1,num_sessions,isesh)
   imagesc(zscore_mtx(tuning_curve_matrix_ref(:,:,isesh)')')
   caxis([-1.5 1.5])
end


