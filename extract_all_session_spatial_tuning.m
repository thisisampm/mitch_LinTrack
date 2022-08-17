function [tuning_curve_mtx_all_no_preference, tuning_curve_mtx_all_left, tuning_curve_mtx_all_right] =...
    extract_all_session_spatial_tuning(session_files, cell_regist_mtx, trace_type)

nbins = 40;
ncells = size(cell_regist_mtx,1); % Each cell reg row represents a tracked cell
nsessions = size(cell_regist_mtx,2); % Each cell reg column represents a session

% Preallocate a cell x bin x session matrix outputs
tuning_curve_mtx_all_no_preference = nan(ncells,nbins,nsessions);
tuning_curve_mtx_all_right = nan(ncells, nbins, nsessions); 
tuning_curve_mtx_all_left = nan(ncells, nbins, nsessions);


for isession = 1:nsessions
    load(session_files{isession},'behavior_mtx','traces'); % Load the current session data
    % Select appropriate 2d trace matrix from the 3 choices.
    switch trace_type
        case 1
            trace_mtx = traces(:,:,1);
        case 2
            trace_mtx = traces(:,:,2);
        case 3
            sampling_rate = 20;
            z_threshold = 2;
            trace_mtx = extract_binary(traces(:,:,3),sampling_rate,z_threshold);
    end
    running       = behavior_mtx(:,5) > 5; % Speed threshold over 5 cm/s
    running_right = behavior_mtx(:,5) > 5 & isolate_direction(behavior_mtx(:,2),'right'); % Speed and direction threshold
    running_left  = behavior_mtx(:,5) > 5 & isolate_direction(behavior_mtx(:,2),'left');
    
    % Calculate the current session spatial tuning given restrictions
    [tuning_curve_no_preference] = extract_1D_spatial_tuning(behavior_mtx(running,:),trace_mtx(:,running),trace_type,nbins);
    [tuning_curve_right, ~, ~] = extract_1D_spatial_tuning(behavior_mtx(running_right,:), trace_mtx(:,running_right), trace_type, nbins);
    [tuning_curve_left, ~, ~]  = extract_1D_spatial_tuning(behavior_mtx(running_left,:), trace_mtx(:,running_left), trace_type, nbins); 

    % Determine the rows of the cell reg matrix that correspond to the
    % current session cells
    active_cells = cell_regist_mtx(:,isession) > 0;
    active_cell_order = cell_regist_mtx(active_cells,isession);
    
    % Assign current session tuning curves to their appropriate position in all session matrx
    tuning_curve_mtx_all_no_preference(active_cells,:,isession) = tuning_curve_no_preference(active_cell_order,:);
    tuning_curve_mtx_all_left(active_cells,:,isession) = tuning_curve_left(active_cell_order,:); 
    tuning_curve_mtx_all_right(active_cells,:,isession) = tuning_curve_right(active_cell_order,:);
end
    
    
    
