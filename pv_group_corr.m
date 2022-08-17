%% Function that takes group and calculates the population vectors correlations for each mouse

function [pv_group_mtx] = pv_group_corr(mouse_group_folder, varargin)

if isempty(varargin)
    reference_session = 1;
else
    reference_session = varargin{1};
end

mouse_files = get_folder_paths_all(mouse_group_folder,1); % Get the folder for each individual mouse in a group
G7_files = contains(mouse_files, 'G7'); % Index to exclude data from an viral GCaMP7f mice
mouse_files = mouse_files(~G7_files); % Exclude GCamP7f mice
n_session = numel(get_file_paths_all(mouse_files{1}))-1; % Calculate the number of sessions as the number of files in a mouse folder minus one CellReg file
n_mice = numel(mouse_files); % Calculate number of mice

pv_group_mtx = nan(n_session, n_mice); % Preallocate matrix to store average PV correlations with each mouse's data in a column

% For each mouse calculate the average PV correlation between session 1 and
% all other sessions 
for i = 1:n_mice
    pv_group_mtx(:,i) = pv_session_corr(mouse_files{i}, reference_session);
end
    
    
    
    