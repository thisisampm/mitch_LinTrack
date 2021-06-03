function dot_plot(behavior_mtx, traces, neuron)
% plots of activity over space, only from included trials

% figure
figure; 

%% overlay
subplot(1,2,1); hold on; 

% trials only
trl_idx = ~isnan(behavior_mtx(:,4));

% neuron index
n_idx = (traces(neuron,:,1)==1)';

% plot behavior
plot(behavior_mtx(trl_idx,2), behavior_mtx(trl_idx,3), '-', 'color', .7.*[1 1 1]); 

% plot spike locations
plot(behavior_mtx(trl_idx & n_idx, 2), behavior_mtx(trl_idx & n_idx, 3), 'r.'); 

% aesthetics
axis equal; maze_outline
ylabel('y position')
xlabel('x position')
legend({'Location', 'Spike'})


%% time

subplot(1,2,2); hold on; 

% trials only
unq_trials = unique(behavior_mtx(~isnan(behavior_mtx(:,4)),4))';

%iterate through trials
for itrl = unq_trials

    % trial index
    local_trl_idx = behavior_mtx(:,4)==itrl;
    
    % plot behavior
    plot(behavior_mtx(local_trl_idx,2), behavior_mtx(local_trl_idx,1), '-', 'color', .7.*[1 1 1]); 

    % plot spike locations
    plot(behavior_mtx(local_trl_idx & n_idx, 2), behavior_mtx(local_trl_idx & n_idx,1), 'r.'); 

end

% aesthetics
set(gca, 'TickLength', [0 0])
fulltime = [min(behavior_mtx(trl_idx,1)) max(behavior_mtx(trl_idx,1))];
ylim([fulltime(1)-diff(fulltime)*.1 fulltime(2)+diff(fulltime)*.1])
xlim([-.1 1.1]);
ylabel('Time (s)')
xlabel('x position')



%%

% figure; position
set(gcf, 'Position', [168         714        1235         345])
