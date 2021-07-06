function mitch_trial_finder(behavior_mtx)
figure; hold on
plot3(behavior_mtx(:,2), behavior_mtx(:,1), behavior_mtx(:,3), '-', 'color', 0.7.*[1 1 1])
for itrl = unique(behavior_mtx(~isnan(behavior_mtx(:,4)),4))'
    trl_x = behavior_mtx(behavior_mtx(:,4)==itrl,2);
    trl_t = behavior_mtx(behavior_mtx(:,4)==itrl,1);
    trl_y = behavior_mtx(behavior_mtx(:,4)==itrl,3);
    
    % excluded trials plotted with black
    %if ismember(itrl,excluded_trials)
     %   plot3(trl_x, trl_t, trl_y, 'k-', 'linewidth', 2)
   % else
    %    
    %end
    plot3(trl_x, trl_t, trl_y, 'linewidth', 2)
end
ylim([-10 inf])
ylabel('time (s)')
xlabel('x position')
set(gca,'TickLength',[0, 0]); box off;