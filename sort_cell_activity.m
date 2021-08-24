function [sorted_rate_mtx, rate_mtx, rate_mtx_3d] = sort_cell_activity(behavior_mtx, trace_mtx, num_bins)
% plot rate matrix of spatial trial activity for every neuron in traces
% sort cells by location of FR peak
%
% trace mtx is (n,n,1)
%

% minimum number of trials that must have activity for cell to count as
% active
act_trl_min = 3;

% left and right trial index
LRidx = left_right_trial_idx(behavior_mtx);
LRidx = LRidx(~isnan(LRidx));

% unique trials
unq_trials = unique(behavior_mtx(~isnan(behavior_mtx(:,4)),4));
rightward_trials = unq_trials(LRidx==1);
leftward_trials = unq_trials(LRidx==0);
min_LR_trials = min([length(rightward_trials) length(leftward_trials)]);

% compute activity on every trial for every cell
% CHANGE HERE TO DETERMINE IF WANT TO BREAK TRACK UP INTO LR
rate_mtx_3d = sesh_bin_trl_act(behavior_mtx, trace_mtx, num_bins, unq_trials); % LR combined
%rate_mtx_3d = [sesh_bin_trl_act(behavior_mtx, trace_mtx, num_bins, rightward_trials(1:min_LR_trials)) sesh_bin_trl_act(behavior_mtx, trace_mtx, num_bins, leftward_trials(1:min_LR_trials))]; % LR separate


% drop cells that aren't active on at least 3 trials
for ic = 1:size(rate_mtx_3d,3)
   cell_activity = rate_mtx_3d(:,:,ic)~=0;
   num_trl_act = max(sum(cell_activity));
   if num_trl_act < act_trl_min
       rate_mtx_3d(:,:,ic) = nan;
   end
end

% average over trials
rate_mtx = squeeze(nanmean(rate_mtx_3d))';

% sort rows by peak
sorted_rate_mtx = sort_rows_by_peak(rate_mtx);

% plot sorted (normed)
figure; imagesc(norm_mtx(sorted_rate_mtx')')







