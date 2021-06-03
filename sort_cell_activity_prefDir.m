function [sorted_rate_mtx, rate_mtx] = sort_cell_activity_prefDir(behavior_mtx, trace_mtx, num_bins)
% plot rate matrix of spatial trial activity for every neuron in traces
% sort cells by location of FR peak
%
% trace mtx is (n,n,1)
%

% unique trials
unq_trials = unique(behavior_mtx(~isnan(behavior_mtx(:,4)),4));

% travel direction
dir_vect = trial_direction(behavior_mtx);

% compute activity on every trial for every cell
rate_mtx_3d = sesh_bin_trl_act(behavior_mtx, trace_mtx, num_bins, unq_trials);

% find preferred direction for each cell
[~, mean_corrs] = ALL_pairwise_trial_corrs(behavior_mtx, rate_mtx_3d);
[~,pref_dir] = max(mean_corrs(:,2:3),[],2); pref_dir = pref_dir-1;

% average over preferred direction trials
rate_mtx = nan(length(unq_trials),num_bins);
for ic = 1:size(rate_mtx_3d,3)
    rate_mtx(ic,:) = mean(rate_mtx_3d(dir_vect==pref_dir(ic),:,ic));
end

% sort rows by peak
sorted_rate_mtx = sort_rows_by_peak(rate_mtx);

% plot sorted (normed)
figure; imagesc(norm_mtx(sorted_rate_mtx')')







