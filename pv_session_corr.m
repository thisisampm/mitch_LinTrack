%% Function that takes a mouse folder, loads data, generates and compares

function [pv_corr_mtx] = pv_session_corr(mouse_folder, varargin)

% If no reference session specified with varargin set the reference session
% to first session
if isempty(varargin)
    reference_session = 1;
else
    reference_session = varargin{1};
end

%% Collect the files to analyze from the mouse folder
session_files = get_file_paths_all(mouse_folder); % Get all the session files in the folder
cell_reg_index = contains(session_files, {'cellReg', 'cell_regist'},'IgnoreCase',1); % Index the files for session data vs cellReg file
cell_reg_file = session_files(cell_reg_index);
session_files = session_files(~cell_reg_index);

%% Load and calculate the mean PV each spatial bin during each session

pv_cell = cell(size(session_files)); % Preallocate a cell to hold the average population vector for each spatial bin
for i = 1:numel(session_files)
    load(session_files{i}) % Load the appropriate file containing behaviour and Ca trace data
    [~, pv_cell{i}] = pv_trial(behavior_mtx, traces); % Calculate and store mean PV
end

%% Load cellReg cell idenity matrix
variableInfo = who('-file',cell_reg_file{1}); % Check the cell reg file for variables
        if ismember('cell_regist_mtx', variableInfo)
            load(cell_reg_file{1}, 'cell_regist_mtx') % The name saved by Adam maybe instead of the whole cell reg file??
            cell_reg_mtx = cell_regist_mtx;
        elseif ismember('cell_registered_struct', variableInfo)
            load(cell_reg_file{1}, 'cell_registered_struct') % The name of the Cell Reg output structure
            cell_reg_mtx = cell_registered_struct.cell_to_index_map;
        else
            error('bad cell_regist.mat file')
        end


%% Get the common session between sessions and compare their population vectors

pv_corr_mtx = nan(size(session_files)); % Preallocate matrix to store the mean PV correlation between reference session and each other session.

for i = 1:numel(session_files)
    common_cells = cell_reg_mtx(:,reference_session) > 0 & cell_reg_mtx(:,i) > 0; % Logical index of cells registered as the same same active on reference and comparison session
    ref_cell_index = cell_reg_mtx(common_cells, reference_session); % Ordered index of active cells in reference session
    comparison_cell_index = cell_reg_mtx(common_cells, i); % Ordered index of active cells in comparison session
    pv_ref = pv_cell{reference_session}(ref_cell_index,:); % Population vector activity on reference session for all cells common to comparison session
    pv_comparison = pv_cell{i}(comparison_cell_index,:); % Population vector activity on comparison session for all cells common to comparison session
    pv_correlations = corr(pv_ref,pv_comparison); % Calculate the correlations between the population vectors for each of the spatial bins between reference and comparison
    pv_corr_mtx(i) = mean(diag(pv_correlations),'omitnan'); % Comparisons between same spatial bin lie along the diagonal of the corr mtx.
end


