function rate_mtx_3d = sesh_bin_trl_act(behavior_mtx, trace_mtx, num_bins, trial_nums)
% outputs four dimensional activity mtx (trial, bin, neuron)
%
% trace mtx is (n,n,1)

% preallocate rate mtx
rate_mtx_3d = nan(length(trial_nums), num_bins, size(trace_mtx,1));

% compute average trial rate vect for every neuron
for ic = 1:size(trace_mtx,1)
    
    % load average activity
    rate_mtx_3d(:,:,ic) = binned_trl_activity(behavior_mtx, trace_mtx(ic,:)', num_bins, trial_nums);
    
end