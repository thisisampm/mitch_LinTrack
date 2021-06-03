function [fpaths, unq_subj] = get_file_paths_targeted_III(folderpath, session_strings, subj_strings, exclusion_strings)
% returns files on 'folderpath' that contain all the strings in 
% 'session_strings', but ONLY from subjects who have a least 1 file path 
% containing each of the strings in 'subj_strings'. Excludes all paths
% containing strings in 'exclusion_strings'.
%
% e.g., path_comp_II('C:\path', {'probe-02'}, {'newlearn05','probe_01d'})
%

% all paths under umbrella
fpaths = get_file_paths_all(folderpath);

% constrain fpaths by subject strings
if ~isempty(subj_strings)
    
    % get all subjects for each of the subject strings
    unq_subj_hold = cell(length(subj_strings),1);
    for isubstr = 1:length(subj_strings)
        
        fpath_match = fpaths(contains(fpaths, subj_strings{isubstr}));
        
        if isempty(fpath_match)
            disp('no paths')
            fpaths = [];
            return
        end
        
        unq_subj_hold{isubstr} = unique_subjects(fpath_match);
    end
    unq_subj = mintersect(unq_subj_hold, 'rows');

    % constrain fpaths by unique subject list
    idx_hold = zeros(length(fpaths), length(unq_subj));
    for isubj = 1:size(unq_subj,1)
        idx_hold(contains(fpaths, unq_subj(isubj,:)),isubj) = 1;
    end
    fpaths = fpaths(sum(idx_hold,2)>0);
    
    if isempty(fpaths)
        disp('no paths')
        fpaths = [];
        return
    end
    
end


% constrain via session strings
for istr = 1:length(session_strings)
    fpaths = fpaths(contains(fpaths, session_strings{istr}));
    
    if isempty(fpaths)
        disp('no paths')
        fpaths = [];
        return
    end
    
end


% constrain via exclusion strings
for istr = 1:length(exclusion_strings)
    fpaths = fpaths(~contains(fpaths, exclusion_strings{istr}));
end





% internal function
function us_out = unique_subjects(file_paths)
% get list of unique subjects

    % find bounds of subject folder
    slash_pos = strfind(file_paths{1},'\');
    slash_pos = [slash_pos(end-1) slash_pos(end)];

    % snip subject id from each path
    us_out = cellfun(@(v) v(slash_pos(1)+1 : slash_pos(2)-1),file_paths, 'UniformOutput', 0);

    % char array
    us_out = cell2mat(us_out);
    
    % unique subjects
    us_out = unique(us_out, 'rows');


