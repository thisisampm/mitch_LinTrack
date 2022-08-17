function [all_tcm_no_direction, all_tcm_left, all_tcm_right, pc_mtx_cell, combined_pc_table, mice, pc_tcm_cell]...
    = extract_1D_spatial_tuning_multi(group_folder)
% Function to create the tuning curve matrices for each mouse in a group

trace_type = 3;
[folder_paths, mice] = get_folder_paths_all(group_folder); % Get all of the individual mice contained within a group folder
% Preallocate cell arrays that will hold the matrix/array produced from
% each mouse
pc_mtx_cell = cell(size(mice));
all_tcm_no_direction = cell(size(mice));
all_tcm_left = cell(size(mice));
all_tcm_right = cell(size(mice));
pc_tcm_cell = cell(size(mice));

% For each mouse load its tuning curves and place cell matrix
for imouse = 1:numel(mice)
    file_paths = get_file_paths_all(folder_paths{imouse}); % Cell individual files within mouse folder
    reg_id = contains(file_paths,{'cell_reg','cellreg'},'IgnoreCase',1); % identify the cell reg file
    session_files = file_paths(~reg_id); % Session files containing behaviour and calcium data
    reg_file = file_paths(reg_id); % The cell reg alignment across the session files
    load(reg_file{:},'cell_registered_struct');
    cell_reg_mtx = cell_registered_struct.cell_to_index_map; % Each row contains information about which cells have been aligned across days
    [imouse_pc_mtx, imouse_pc_table] = all_session_place_cell_info(folder_paths{imouse},cell_reg_mtx); % Calculate matrix indicating if a cell has a significant place field in each session
    if ~exist('combined_pc_table','var')
        combined_pc_table = imouse_pc_table;
    elseif exist('combined_pc_table','var')
        combined_pc_table = [combined_pc_table; imouse_pc_table];
    end
    % Calculate tuning curves without and with each directional preference
    [imouse_all_tcm_no_direction, imouse_all_tcm_left, imouse_all_tcm_right] = extract_all_session_spatial_tuning(session_files,cell_reg_mtx,trace_type);
    left_place_field = any(imouse_pc_mtx == 1 | imouse_pc_mtx == 3,2); % Any cell with significant left place field
    right_place_field = any(imouse_pc_mtx == 2 | imouse_pc_mtx == 3,2); % 
    % Try setting all non-significant place fields back to NaNs
    imouse_all_tcm_left(~left_place_field,:,:) = nan;
    imouse_all_tcm_right(~right_place_field,:,:) = nan;
    % Make a combined matrix of place cell tuning in preferred direction
    pc_tcm_imouse = nan(size(imouse_all_tcm_no_direction));
    pc_tcm_imouse(left_place_field,:,:) = imouse_all_tcm_left(left_place_field,:,:);
    pc_tcm_imouse(right_place_field,:,:) = imouse_all_tcm_right(right_place_field,:,:);
    pc_tcm_cell{imouse} = pc_tcm_imouse;
    % Assign outputs into their cell array
    pc_mtx_cell{imouse} = imouse_pc_mtx;
    all_tcm_no_direction{imouse} = imouse_all_tcm_no_direction;
    all_tcm_left{imouse} = imouse_all_tcm_left;
    all_tcm_right{imouse} = imouse_all_tcm_right;
end
    
    
    
    
