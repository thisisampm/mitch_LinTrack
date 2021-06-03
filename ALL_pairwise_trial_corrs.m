function [max_mean_corrs, mean_corrs] = ALL_pairwise_trial_corrs(behavior_mtx, rate_mtx_3d)
% iterates through neurons trial corr matrices to find highest r val for
% each cell
%
% input
% [rate_mtx_3d] = sesh_bin_trl_act(behavior_mtx, trace_mtx, num_bins);
%

% preallocate
mean_corrs = nan(size(rate_mtx_3d,3),3);

% iterate through all neurons
for ic = 1:size(rate_mtx_3d,3)
    
    % compute correlations from all, left, and right trials
    mean_corrs(ic,:) = pairwise_trial_corrs(behavior_mtx, rate_mtx_3d(:,:,ic));
    
end

% maximums
max_mean_corrs = max(mean_corrs(:,2:3),[],2);

