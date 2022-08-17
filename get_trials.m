%% Function to extract the trials information LinTrack behaviour mtx
% MDS 25 Mat 2022

% Input: behavior_mtx
% Outputs: 
% 1. the index of unique trial numbers
% 2. trialx2 matrix with the start and end frame for each trial
% 3. A vector that indicates trial direction 0 = rightward 1 = leftward

function [unique_trials, trial_markers, direction] = get_trials(behavior_mtx)

% Index of the number of trials
unique_trials = unique(behavior_mtx(~isnan(behavior_mtx(:,4)),4));

% Get the start and end frame of each trial and create a vector that
% indicates the direction of each trial
trial_markers = nan(length(unique_trials),2); % Preallocate
direction = nan(length(unique_trials),1); 
for i = 1:numel(unique_trials)
    trial_markers(i,1) = find(behavior_mtx(:,4) == unique_trials(i),1,'first');
    trial_markers(i,2) = find(behavior_mtx(:,4) == unique_trials(i),1,'last');
    % Create direction vector 0 = left-to-right (rightward), 1 =
    % right-to-left (leftward)
    if behavior_mtx(trial_markers(i,1),2) < behavior_mtx(trial_markers(i,2),2)
        direction(i) = 0;
    elseif behavior_mtx(trial_markers(i,1),2) > behavior_mtx(trial_markers(i,2),2)
        direction(i) = 1;
    end    
end

    

