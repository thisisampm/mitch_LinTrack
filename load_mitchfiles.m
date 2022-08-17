function load_mitchfiles(folderpath, varargin)
% finds data located in folder_path and saves it in the mitch_linTrack 
% matlab folder according to mouse and day
%
% folderpath = 'D:\Projects\InProgress\Mitch_LinTrack\data';

destination_path = 'C:\Users\mitch\OneDrive - University of Toronto\MATLAB\mitch_LinTrack\reload';

if nargin == 2
    specific_mouse = varargin{1};
elseif nargin == 3
    specific_mouse = varargin{1};
    dest_group = varargin{2};
end

% group folders
[group_folder_paths, group_folder_names] = get_folder_paths_all(folderpath,1);

% ensure destination data folder exists 
dest_data_fp = [destination_path];
if ~exist(dest_data_fp, 'dir')
    mkdir(dest_data_fp)
end


% iterate through group folders
for igfolder = 1:size(group_folder_paths,1)
    
    % ensure destination group folder exists 
    dest_group_fp = [destination_path '\' group_folder_names{igfolder}];
    if ~exist(dest_group_fp, 'dir')
        mkdir(dest_group_fp)
    end
    
    % get mouse folders within group folder
    [mouse_folder_paths, mouse_folder_names] = get_folder_paths_all(group_folder_paths{igfolder},1);


    % iterate through mouse folders
    for imfolder = 1:size(mouse_folder_paths,1)
        
        % If varagin a specfic mouse, skip the others.
        if exist('specific_mouse','var')
            if ~strcmp(specific_mouse, mouse_folder_names{imfolder})
                continue
            end
        end
        
        % ensure distination mouse folder exists 
        dest_mouse_fp = [destination_path '\' group_folder_names{igfolder} '\' mouse_folder_names{imfolder}];
        if ~exist(dest_mouse_fp, 'dir')
            mkdir(dest_mouse_fp)
        end
        
        
        % get session folders within mouse folder
        [sesh_file_paths, sesh_file_names] = get_folder_paths_all(mouse_folder_paths{imfolder},1);
        
        
        % exceptions
        sesh_file_paths = sesh_file_paths(~contains(sesh_file_names, 'ootprints'));
        sesh_file_names = sesh_file_names(~contains(sesh_file_names, 'ootprints'));
        
        % iterate through origin session folders
        for isfolder = 1:size(sesh_file_paths,1)
            
            % only save session if not already saved
            mat_filename = sesh_file_names{isfolder};
            save_fn = [dest_mouse_fp '\' mat_filename];
            exist([save_fn '.mat'], 'file')
            if ~exist([save_fn '.mat'], 'file')
                
                % load file
                [behavior_mtx, traces] = load_mitchdata(sesh_file_paths{isfolder});

                % save file
                save(save_fn, 'behavior_mtx', 'traces')
                
                % show figures
                drawnow
            end
            
            

        end
    
    end
end