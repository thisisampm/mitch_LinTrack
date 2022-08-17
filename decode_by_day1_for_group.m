function [mean_error, median_error, true_positions, decoded_positions, mouse_ids] = decode_by_day1_for_group(group_folder)
    
%% Get the mouse info
[mouse_folders, mouse_ids] = get_folder_paths_all(group_folder);

n_mice = numel(mouse_folders);
n_sessions = numel(get_file_paths_all(mouse_folders{1}))-1;

%% Preallocate outputs
mean_error = nan(n_sessions-1, n_mice); % Co
median_error = nan(n_sessions-1, n_mice);
true_positions = cell(n_mice,1);
decoded_positions = cell(n_mice,1);

%% Calculate for each mouse in group
for i_mouse = 1:numel(mouse_ids)
    [mean_error(:,i_mouse), median_error(:,i_mouse), true_positions{i_mouse}, decoded_positions{i_mouse}] = ...
        decode_by_day1(mouse_folders{i_mouse});
end
