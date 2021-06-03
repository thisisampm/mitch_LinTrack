function [dir_vect] = trial_direction(behavior_mtx)
% outputs vector of zeros (leftward) and ones (rightward) describing the
% direction of the mouse on that trial

% unique trials
unq_trials = unique(behavior_mtx(~isnan(behavior_mtx(:,4)),4))';

% preallocate
dir_vect = nan(length(unq_trials),1);

% iterate through trials
for itrl = 1:length(unq_trials)
   trial_xpos = behavior_mtx(behavior_mtx(:,4)==unq_trials(itrl),2);
   if trial_xpos(1)<trial_xpos(end)
       dir_vect(itrl) = 0;
   elseif trial_xpos(1)>trial_xpos(end)
       dir_vect(itrl) = 1;
   end
end
