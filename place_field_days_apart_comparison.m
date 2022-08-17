function [] = place_field_days_apart_comparison(pc_mtx)

% Input: pc_mtx where rows represent individual cells and columns session
% enteries represent place cell classification for that cell on that day
% NaN = Not active, 0 = no place field, 1 = left place field , 2 = right
% place field, 3 = place fields in both directions.

n_cells = size(pc_mtx,1);
n_sessions = size(pc_mtx,2);

% Identify all possible n+1,n+3,n+5...n+max_separation comparisons that can be made
sessions = 1:n_sessions; % Row vector indicating session #
session_separation_mtx = sessions - sessions'; % Matrix indicating how far any given session is from any other session
comp_distances = unique(session_separation_mtx(session_separation_mtx > 0)); % Vector indicating all the unique session separation possibilities



% For each comparison distance 
for idist = comp_distances'
    [s1, s2] = find(session_separation_mtx == idist); % Find which sessions correspond to sessions i days apart
    % For each session n days apart, determine how many cells with place
    % fields maintain a place field
    for icomp = 1:length(s1)
        s1_active_cells = ~isnan(pc_mtx(:,s1(icomp)));
        s1_place_cells = pc_mtx(s1_active_cells,s1(icomp)) ~= 0; % Cells in the first session that have a significant place field
        s2_place_cells = pc_mtx(s1_active_cells,s2(icomp)) ~= 0; % Place cell classification in the 2nd session of cells in the 1st session
        cell_overlap = []; % The fraction of place cells that remain place cells 
        
    end
    
end

