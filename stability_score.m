function [stability_score] = stability_score(behavior_mtx, traces)
% Function to compare the stability of place fields in the first half of
% trials to the second half of trials

% Inputs: 
% behavior_mtx: full 6 column behaviour matrix with time, xpos, ypos,
% trial indicator, velocity, and acceleration
% trace_mtx: neuron x timestep x trace type, 2d matrix with Ca activity data matching
% behavior_mtx

% Output:
% SC: 

%%
% Get trial information
[trials, ~, direction] = get_trials(behavior_mtx);

% Cell activity for each cell on each trial
trace_type = 3; % Binarized C_raw
nbins = 40; % Number of spatial bins
rate_mtx_3d = binned_trl_activity_MDS(behavior_mtx, traces, trace_type, nbins, trials); % Trial x bin x neuron Ca activity matrix

% Logical vectors indicating which direction mouse is running a trial
rightward = direction == 0;
leftward =  direction == 1;

% Logical vector indicating which direction a neuron is more often active
active_right = sum((sum(rate_mtx_3d(rightward,:,:) > 0,2))>0); % Returns 1x1xneuron matrix with the number of trials which the cell is active
active_left = sum((sum(rate_mtx_3d(leftward,:,:) > 0,2))>0); % Repeat for left
fraction_active_right = squeeze(active_right)/sum(rightward); % Expressed as a fraction of trials
fraction_active_left = squeeze(active_left)/sum(leftward); % Repeat left
pref_right = fraction_active_right > fraction_active_left; % Logical index, marks cells more frequently active on rightward trials
pref_left = ~pref_right; 

% For each neuron, calculate a stability score as the Pearson correlation 
% of its tuning curve on first half trials compared to its second half of
% trials in its preferred direction
right_trials = find(rightward);
first_half_right = right_trials(1:floor(numel(right_trials)/2));
second_half_right = right_trials(~ismember(right_trials,first_half_right));

left_trials = find(leftward);
first_half_left = left_trials(1:floor(numel(left_trials)/2));
second_half_left = left_trials(~ismember(left_trials,first_half_left));

% Calculate tuning curve as the average activity in a spatial bin
% Calculate tuning for both directions for each neuron
tuning_curves_rightward_first_half = squeeze(mean(rate_mtx_3d(first_half_right,:,:),1,'omitnan'));
tuning_curves_rightward_second_half = squeeze(mean(rate_mtx_3d(second_half_right,:,:),1,'omitnan'));

tuning_curves_leftward_first_half = squeeze(mean(rate_mtx_3d(first_half_left,:,:),1,'omitnan'));
tuning_curves_leftward_second_half = squeeze(mean(rate_mtx_3d(second_half_left,:,:),1,'omitnan'));

% Extract the tuning curve for the preferred direction for each neuron
tuning_curve_pref_dir_first_half = nan(size(tuning_curves_leftward_first_half));
tuning_curve_pref_dir_first_half(:,pref_right) = tuning_curves_rightward_first_half(:,pref_right);
tuning_curve_pref_dir_first_half(:,pref_left) = tuning_curves_leftward_first_half(:,pref_left);

tuning_curve_pref_dir_second_half = nan(size(tuning_curves_leftward_second_half));
tuning_curve_pref_dir_second_half(:,pref_right) = tuning_curves_rightward_second_half(:,pref_right);
tuning_curve_pref_dir_second_half(:,pref_left) = tuning_curves_leftward_second_half(:,pref_left);

%% 
% Pearson correlation of first half tuning curves to second half
correlation_matrix = corr(tuning_curve_pref_dir_first_half,tuning_curve_pref_dir_second_half);
% Extract correlations between tuning curves of the same cell (diagonal of
% matrix returned by corr) and Fisher transform to give stability score
corr_first_to_second_half = diag(correlation_matrix);
nancorrs = find(isnan(corr_first_to_second_half)); % Identify which cells have nan correlations
% Check if they lack activity on either half of the trials
first_half_no_activity = all(tuning_curve_pref_dir_first_half(:,nancorrs) == 0);
second_half_no_activity = all(tuning_curve_pref_dir_second_half(:,nancorrs) == 0);
% If no activity on either half of trials, set correlation to 0
% (uncorrelated)
corr_first_to_second_half(nancorrs(first_half_no_activity | second_half_no_activity)) = 0;

% Fisher transform the correlations
stability_score = atanh(corr_first_to_second_half); 



