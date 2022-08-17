%% Function that identifies which cells are place cells across days.

function pc_map = place_cell_to_index_map(mouse_folder)

%% Get mouse cell reg file and session files
files = get_file_paths_all(mouse_folder);
reg_idx = contains(files, {'cell_reg','cellReg'}, "IgnoreCase", 1);
reg_file = files{reg_idx};
session_files = files(~reg_idx);

load(reg_file); % Load the cell reg file

% Most folders contain the whole cell reg output structure, some are
% already just the cell_to_index_map so need to appropriately extract map.
if ~exist('cell_regist_mtx', 'var')
    cell_regist_mtx = cell_registered_struct.cell_to_index_map;
end

%% Map the session place cell index to their corresponding position in output matrix

pc_map = nan(size(cell_regist_mtx)); % Preallocate an cell x session matrix (same format as cell reg)

for i = 1:numel(session_files)
    load(session_files{i},'place_cell_mtx'); % Load place cell data
    active_cells = cell_regist_mtx(:,i) > 0; % Identify cell's active in ith session
    active_cell_order = cell_regist_mtx(active_cells,i); % 
    place_cell_index = place_cell_mtx(:,2) == 1; % Logical index of place cell positions in current session cnfme order
    pc_map(active_cells,i) = place_cell_index(active_cell_order); % Place indicator at output position if a cell meets place cell criteria on a given day
end