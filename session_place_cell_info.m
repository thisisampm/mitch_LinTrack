function [session_pc_mtx] = session_place_cell_info(behavior_mtx,traces)
% function to determine the place tuning for each cell in each session
% based on the direction of running

%% Initialize
running_ts = behavior_mtx(:,5) > 5; % Index out only when the mouse is running faster than 5 cm/s
left = isolate_direction(behavior_mtx(:,2),'left');
right = isolate_direction(behavior_mtx(:,2),'right');

%% Calculate spatial information, place field significance
[si_left, pc_idx_left, pc_peak_left, pf_width_left] = information_score(behavior_mtx, traces, running_ts&left, 'left');
[si_right, pc_idx_right, pc_peak_right, pf_width_right] = information_score(behavior_mtx, traces, running_ts&right, 'right');

%% Combined data
ncells = size(traces,1);
% 0 = no place field
% 1 = place field left
% 2 = place field right
% 3 = place field both directions
session_pc_mtx = nan(ncells,1);
session_pc_mtx(~(pc_idx_left | pc_idx_right)) = 0;
session_pc_mtx(pc_idx_left) = 1;
session_pc_mtx(pc_idx_right) = 2;
session_pc_mtx(pc_idx_left & pc_idx_right) = 3;

session_pc_mtx = [session_pc_mtx, si_left, si_right, pc_peak_left, pc_peak_right, pf_width_left, pf_width_right];