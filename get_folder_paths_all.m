function [item_paths, folder_names] = get_folder_paths_all(folderpath, level)
%iterates through each folder in path and returns path incl folder

item_paths = cell(1);

%items in path
file_list = dir(folderpath);
for iitem = 1:length(file_list)
    current_sesh = file_list(iitem).name;
    
    %omited folder and file names (folders named 'old' are invisible)
    if strcmp(current_sesh, '.') || strcmp(current_sesh, '..') || strcmp(current_sesh, 'old')
        continue
    end

    if isfolder([folderpath '\' file_list(iitem).name])
        item_paths = [item_paths; {[folderpath '\' file_list(iitem).name]}];
        item_paths_hold = get_folder_paths_all([folderpath '\' file_list(iitem).name]);
        item_paths = [item_paths; item_paths_hold];
    end

end

%remove empty cells
item_paths = item_paths(find(~cellfun(@isempty, item_paths)));



% only include folders from specified level
if exist('level', 'var')
    
    % num levels in input folder
    origin_level = length(strfind(folderpath , '\'));
    
    % num levels in each identified folder
    pathlevels = nan(size(item_paths,1),1);
    for ipath = 1:size(item_paths,1)
        pathlevels(ipath) = length(strfind(item_paths{ipath} , '\'));
    end

    % only keep paths with the specified level
    item_paths = item_paths(pathlevels == origin_level+level);
    
end


% get folder names
folder_names = cell(size(item_paths,1),1);
for ipath = 1:size(item_paths,1)
    slash_pos = strfind(item_paths{ipath} , '\');
    folder_names{ipath} = item_paths{ipath}(slash_pos(end)+1:end);
end

