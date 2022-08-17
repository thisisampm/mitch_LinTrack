% Function to update the velocity vector

function update_velocity(linear_track_folder)

mouse_folders = get_folder_paths_all(linear_track_folder,2); % Get all the individual mouse folders below the group folders
track_length = 91; % Track is 91 cm long

for imouse = 1:numel(mouse_folders)
    file_names = get_file_paths_all(mouse_folders{imouse}); % Get session files
    cell_reg_file_index = contains(file_names, {'cellReg','Cell_Reg'},'IgnoreCase',1); % Find the cell reg file
    file_names = file_names(~cell_reg_file_index); % Exclude cell reg file
    for isesh = 1:numel(file_names)
        clearvars behavior_mtx
        disp(['Updating velocity vector for ', file_names{isesh}])
        load(file_names{isesh},'behavior_mtx','traces','place_cell_mtx') % Load current behavior_mtx file
        velocity = nan(length(behavior_mtx),1); % Preallocate velocity vector
        velocity(1) = 0; % Set time(0) velocity to 0
        x2 = behavior_mtx(2:end,2); x1 = behavior_mtx(1:end-1,2); % Position data from behavior_mtx
        y2 = behavior_mtx(2:end,3); y1 = behavior_mtx(1:end-1,3);
        t2 = behavior_mtx(2:end,1); t1 = behavior_mtx(1:end-1,1); % Time data from behavior_mtx
        % Calculate the velocity as the change in position over change in time
        % The position data is recorded as the relative distance along the
        % track so need to multiple by the track length in cm to get cm/s
        velocity(2:end) = track_length*sqrt((x2 - x1).^2 + (y2 - y1).^2)./(t2 - t1); 
        velocity = smooth(velocity,1/mode(diff(t2))); % Smooth with a 1 second window, corresponding to 1/framerate
        behavior_mtx(:,5) = velocity;
        disp(['Saving ', file_names{isesh}]);
        save(file_names{isesh},'behavior_mtx','traces','place_cell_mtx');
    end
end

        