function [sorted_rate_mtx, rate_mtx] = sort_cell_activity(behavior_mtx, trace_mtx, num_bins)
% plot rate matrix of spatial trial activity for every neuron in traces
% sort cells by location of FR peak
%
% trace mtx is (n,n,1)
%

% minimum number of trials that must have activity for cell to count as
% active
act_trl_min = 3;


% unique trials
unq_trials = unique(behavior_mtx(~isnan(behavior_mtx(:,4)),4));

% compute activity on every trial for every cell
rate_mtx_3d = sesh_bin_trl_act(behavior_mtx, trace_mtx, num_bins, unq_trials);

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







