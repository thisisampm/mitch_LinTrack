function LRidx = left_right_trial_idx(behavior_mtx)
% index of leftward and rightward traveling trials
% start on left (rightward) == 0
% start on right (leftward) == 1

% preallocate
LRidx = nan(length(unique(behavior_mtx(~isnan(behavior_mtx(:,4)), 4))),1);



% iterate through trials
for itrl = unique(behavior_mtx(~isnan(behavior_mtx(:,4)), 4))'

    % id position of first bin
    first_xpos = behavior_mtx(behavior_mtx(:,4)==itrl, 2);
    first_xpos = first_xpos(1);
    
    % if less than .5, rightward
    if first_xpos < 0.5
        LRidx(itrl)=0;
    elseif first_xpos > 0.5
        LRidx(itrl)=1;
    end
    
end

