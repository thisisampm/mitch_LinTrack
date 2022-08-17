function [all_session_pc_mtx, pc_table] = all_session_place_cell_info(mouse_folder, varargin)



%% Load data
files = get_file_paths_all(mouse_folder); % Access all the files within the mouse folder
reg_idx = contains(files, {'cell_reg','cellReg'},'IgnoreCase',1); % Identify the cell_reg file
session_files = files(~reg_idx);
reg_file = files(reg_idx);

if ~isempty(varargin)
    cell_reg_mtx = varargin{1};
else
    load(reg_file{:},'cell_registered_struct');
    cell_reg_mtx = cell_registered_struct.cell_to_index_map;
end
%% Get mouse info
names = regexp(mouse_folder, filesep, 'split');
mouse = names{end};
group = names{end - 1};

%% Determine place cell classification for each cell on each active day

% Preallocate a neuron x session mtx to store place cell classification 
all_session_pc_mtx = nan(size(cell_reg_mtx));
% 0 = no place field
% 1 = place field left
% 2 = place field right
% 3 = place field both directions
si_left = nan(size(cell_reg_mtx));
si_right = nan(size(cell_reg_mtx));
peak_left = nan(size(cell_reg_mtx));
peak_right = nan(size(cell_reg_mtx));
width_left = nan(size(cell_reg_mtx));
width_right = nan(size(cell_reg_mtx));
% for each session calculate and store place cell class
fprintf('Calculating place cell statistics for %s mouse %s\n',group,mouse);
for isesh = 1:numel(session_files)
    fprintf('Starting session %i/%i\n',isesh,numel(session_files));
    load(session_files{isesh},'behavior_mtx','traces')
    isesh_pc_mtx = session_place_cell_info(behavior_mtx,traces); % Generate session place cell matrix
    active_cells = cell_reg_mtx(:,isesh) > 0; % Logical index indicating if the cell corresponding to the cell_reg row is active that session
    cell_order = cell_reg_mtx(active_cells,isesh);
    all_session_pc_mtx(active_cells,isesh) = isesh_pc_mtx(cell_order,1);
    si_left(active_cells,isesh) = isesh_pc_mtx(cell_order,2);
    si_right(active_cells,isesh) = isesh_pc_mtx(cell_order,3);
    peak_left(active_cells,isesh) = isesh_pc_mtx(cell_order,4);
    peak_right(active_cells,isesh) = isesh_pc_mtx(cell_order,5);
    width_left(active_cells,isesh) = isesh_pc_mtx(cell_order,6);
    width_right(active_cells,isesh) = isesh_pc_mtx(cell_order,7);
end



%% Create a table

ncells = size(cell_reg_mtx,1);
nsesh = size(cell_reg_mtx,2);

group_vec = categorical(string(repelem(group,ncells*nsesh,1)));
mouse_vec = categorical(string(repelem(mouse,ncells*nsesh,1)));
session_vec = repmat(1:nsesh,ncells,1);
session_vec = session_vec(:);
cell_vec = repmat((1:ncells)',nsesh,1);

pc_table = table(group_vec, mouse_vec, session_vec, cell_vec, all_session_pc_mtx(:),...
    si_left(:),si_right(:),peak_left(:),peak_right(:),width_left(:),width_right(:),...
    'VariableNames',{'Group','Mouse','Session','Cell','Place_cell','si_left','si_right',...
    'peak_left','peak_right','width_left','width_right'});

