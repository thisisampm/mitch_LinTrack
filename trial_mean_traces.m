function trl_means = trial_mean_traces(behavior_mtx, traces, traces_type)
% outputs matrix with mean traces on each trial

unique_trials = unique(behavior_mtx(~isnan(behavior_mtx(:,4)),4))';
trl_means = nan(length(unique_trials), size(traces,1)); 
for itrl = unique_trials 
    trl_idx = behavior_mtx(:,4)==itrl; 
    trl_means(itrl,:) = nanmean(traces(:,trl_idx,traces_type),2);
end