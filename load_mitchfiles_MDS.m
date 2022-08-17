function load_mitchfiles_MDS(folderpath, mouse_id, group_folder)
% finds data located in folder_path and saves it in the mitch_linTrack 
% matlab folder according to mouse and day and spec
%
% folderpath = 'D:\Projects\InProgress\Mitch_LinTrack\data';

% Destination path for all preprocessed linear track data
destination_path = 'C:\Users\mitch\OneDrive - University of Toronto\MATLAB\mitch_LinTrack\reload';

% Check that there is a group folder matching the input group folder
% If not, it's probably a typo and can create the folder outside of this
dest_group_folder = fullfile(destination_path, group_folder);
if ~exist(dest_group_folder, 'dir')
    error('Destination folder does not exist. Check spelling or create folder.')
end

% Check that there is a group subdirectory for the mouse and if not create
% one
mouse_dest_folder = fullfile(dest_group_folder, mouse_id);
if ~exist(mouse_dest_folder, 'dir')
    mkdir(mouse_dest_folder);
    disp(['Created mouse folder: ', mouse_dest_folder])
end

% Get session folder names for this mouse
mouse_fp = fullfile(folderpath, group_folder, mouse_id);
[session_folder_paths, session_folder_names] = get_folder_paths_all(mouse_fp, 1);

%% Copy the cell reg file
footprints_path = session_folder_paths(contains(session_folder_paths, 'footprints','IgnoreCase',1)); % Separate the footprints path to get cell reg file later
reg_source_path = get_file_paths_targeted(footprints_path{1}, 'cellReg');
reg_source_path = reg_source_path{1}; % Pull it out of cell
[~, reg_file_name, reg_ext] = fileparts(reg_source_path);
reg_dest_path = fullfile(mouse_dest_folder, [reg_file_name, reg_ext]);

if ~exist(reg_dest_path,'file')
    disp(['Copied Cell Reg file: ' reg_source_path, ' to ' reg_dest_path]);
    copyfile(reg_source_path, reg_dest_path);
end



%% Copy the trace and behaviour data for each session

% Omit the footprints folder, and other extraneous folders.
session_folder_paths = session_folder_paths(~contains(session_folder_names, {'Exclude','ootprints','Pretraining'}));
session_folder_names = session_folder_names(~contains(session_folder_names, {'Exclude','ootprints','Pretraining'}));

for idir = 1:size(session_folder_names)
    mat_filename = session_folder_names{idir};
    save_fp = fullfile(mouse_dest_folder, mat_filename);
    if ~exist([save_fp '.mat'], 'file')
        disp(['Loading: ', session_folder_paths{idir}]);
        % load file
        [behavior_mtx, traces] = load_mitchdata(session_folder_paths{idir});
        % Save file
        save(save_fp, 'behavior_mtx', 'traces')
        disp(['Saved mat file: ',save_fp]);
        % show figures
        drawnow
    end
end
