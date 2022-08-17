function [tuning_curve_matrix_ref, tuning_curve_mtx_all] = all_cell_tuning(session_files, cell_regist_mtx, reference_session)
% plot the tuning curves of cells that are common across all session files
% sort rows by peak with reference to session number sort_by_sesh

% details
num_spatial_bins = 40;
traces_type = 2; % S, C, RAW

% input
num_sessions = size(session_files,1);

% iterate through sessions computing tuning curves for all cells
tuning_curve_mtx_all = nan(size(cell_regist_mtx,1), num_spatial_bins, num_sessions);
for isesh = 1:num_sessions
    
    % print update
    disp(['Starting session ' num2str(isesh)])
    
    %If session dropped, cell_reg matrix has NaNs, skip to next session
    if isnan(sum(cell_regist_mtx(:,isesh)))
        disp(['Skipping session ', num2str(isesh), 'due to NaN in Cell Registration indicating low N session']); % If a session had low N neurons, I registered the other sessions and inserted a column of NaNs for that session
        continue
    end
    
    % load session
    load(session_files{isesh}, 'behavior_mtx', 'traces')
    
    % load tuning curves
    % inside sort_cell_activity you can set to use L and R runs seperately,
    % may need to change some preallocated matrices slightly
    [~, unsorted_session_tuning_curves] = sort_cell_activity(behavior_mtx, traces(:,:,traces_type), num_spatial_bins);
    session_name = regexp(session_files{isesh},filesep,'split'); % Split the filename into the parts at the file seperator.
    title([session_name{end-2}, ' ',session_name{end-1}, ' ',session_name{end},  ' Spatial Tuning']) % Use the separated file parts to title each session plot
    
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
[mouse_id] = fileparts(session_files{1});
sgtitle([mouse_id ' Spatial Tuning Sorted by Day 1 Activity'])

