function [mean_corr, corr_matrix] = pairwise_trial_corrs(behavior_mtx, trial_rate_mtx)
% finds correlations between every trial pair
% input has trials as rows and columns as spatial bins
% [trial_rate_mtx] = sesh_bin_trl_act(behavior_mtx, trace_mtx, num_bins);
%
% output is the mean correlation and a correlation matrix
%
% mean_corr contains 3 elements: overall, leftward trials only, and 
% rightward trials only

% neurons are columns 
trial_rate_mtx = trial_rate_mtx';

% unique trials
unq_trials = unique(behavior_mtx(~isnan(behavior_mtx(:,4)),4))';

% preallocate
corr_matrix = nan(length(unq_trials));

% iterate through trials
for i1 = 1:length(unq_trials)
    for i2 = 1:length(unq_trials)
        
        if i1==i2
            continue
        end
        
        corr_matrix(i1,i2) = corr(trial_rate_mtx(:, i1), trial_rate_mtx(:, i2));
       
    end
end

% replace nans with zeros
corr_matrix(isnan(corr_matrix)) = 0;

% trial directions
dir_vect = trial_direction(behavior_mtx);

% mean_corr
mean_corr(1) = nanmean(corr_matrix(:));
corm_hold = corr_matrix(dir_vect==0,dir_vect==0);
mean_corr(2) = nanmean(corm_hold(:));
corm_hold = corr_matrix(dir_vect==1,dir_vect==1);
mean_corr(3) = nanmean(corm_hold(:));



