%% Function that takes a linear track mouse folder and outputs quality metrics
% Mitch de Snoo - Mar 2022

function [QC_table] = quality_control(mouse_folder)

%% Get the relevant data
files = get_file_paths_all(mouse_folder); % Get a list of the files in the mouse directory
reg_file_idx = contains(files, {'cellreg' ,'cell_reg'}, 'IgnoreCase', true); % Identify which file corresponds to the CellReg data
session_files = files(~reg_file_idx); % Sessions are all files that are not the CellReg file
reg_file = files{reg_file_idx}; % CelLReg file
[~,mouse_id,~] = fileparts(mouse_folder); % 
% Get the group id for labeling
folder_above = get_folder_paths_all([mouse_folder,filesep,'..']); % 
folder_above = regexp(folder_above{1},filesep,'split');
group_id = folder_above{1};


%% Calculate quality control metrics
% Create a session x quality metric matrix
    % Col 1: Day #
    % Col 2: n total cells detected
    % Col 3: n session cells detected
    % Col 4: n Trials performed
    % Col 5: fraction of total cells active in session
    % Col 6: fraction of cells that overlap with day 1
    % Col 7: n place cells
    
QC = nan(numel(session_files), 7); % Preallocate session x QC metric matrix
load(reg_file); % Load the cell reg file

% Most folders contain the whole cell reg output structure, some are
% already just the cell_to_index_map so need to appropriately extract map.
if ~exist('cell_regist_mtx', 'var')
    cell_regist_mtx = cell_registered_struct.cell_to_index_map;
end

n_total_cells_detected = size(cell_regist_mtx,1); % Number of rows in the Cell reg matrix reflects the total number of cells detected across all sessions

% For each session calculate metrics
for i = 1:numel(session_files)
    load(session_files{i},'behavior_mtx','traces','place_cell_mtx')
    ncells = size(traces,1); % Each row of the traces mtx represents one cell
    ntrials = max(unique(behavior_mtx(~isnan(behavior_mtx(:,4)),4))); % Col 4 of behaviour matrix labels when engaged in a trial
    p_active = ncells/n_total_cells_detected; % Fraction of total that are active in current session
    p_day1_overlap = sum(cell_regist_mtx(:,1) & cell_regist_mtx(:,i))/ncells; % Fraction of i session cells that are also active on session 1
    p_place_cells = sum(place_cell_mtx(:,2))/ncells; % Col 2 in place_cell_mtx indicates if place cell threshold met
    dayn = i*2-1;
    QC(i,:) = [dayn, n_total_cells_detected, ncells, ntrials, p_active, p_day1_overlap, p_place_cells];
end

%% Output as table
QC_table = array2table(QC,'VariableNames',{'Day','Total Cells','Session cells', 'n Trials', '% Active Cells', '% overlap with session 1', '% place cells'});
mouse_id_vector = categorical(cellstr(repelem(mouse_id, size(QC_table,1),1)));
group_id_vector = categorical(cellstr(repelem(group_id, size(QC_table,1),1)));
QC_table = addvars(QC_table,group_id_vector,mouse_id_vector, 'Before',1,'NewVariableNames',{'Group','Mouse'});
return
