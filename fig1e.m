% Create a figure that shows a section of calcium traces aligned to
% behaviour

% Load an example session
load reload\exercise\214-5\2022-03-15.mat
traces(:,:,3) = (norm_mtx(traces(:,:,3)')');

% Determine place cells
inclusion_vector = behavior_mtx(:,5) > 5; % Speed threshold of 5 cm/s
[si, pc_idx] = information_score(behavior_mtx, traces, inclusion_vector, 'none'); % Calculate spatial information and place cell
%%
% Pick the top 15 spatially informative place cells
place_cells = find(pc_idx);
n_cells_to_plot = 20;
[~,top_si_pc] = maxk(si(pc_idx),n_cells_to_plot);
pc_to_plot = place_cells(top_si_pc);
pc_to_plot = sort(pc_to_plot);

% Get trial times
[trials, trial_markers] = get_trials(behavior_mtx);

% Plot the first 5 trials
start_trial = 5;
end_trial = 15;
duration_to_plot = min(trial_markers(start_trial,:)):max(trial_markers(end_trial,:));

trace_type = 3;
[ax,ca_lines] = plot_traces_offset(norm_mtx(traces(pc_to_plot,duration_to_plot,trace_type)')');
%[ax,ca_lines] = plot_trial_patch(behavior_mtx(duration_to_plot,:),(traces(pc_to_plot,duration_to_plot,trace_type)')');


% plot a line indicating 10 seconds
frame_rate = 20;
last_frame = numel(duration_to_plot);
seconds_scale_bar_length = 20; % 20 seconds
hold(ax,'on')
plot(ax,[last_frame - seconds_scale_bar_length*frame_rate, last_frame],[-.1 -.1], 'color','k','LineWidth',2)
set(ca_lines,'color',[0 0 0 0.25])
% Plot a df/F scale bar
plot([last_frame last_frame], [-.1 0.9],'color','k','LineWidth',2);

%Aesthetics
ax.YTick = [0 n_cells_to_plot/2];
ax.YTickLabel = {'1',num2str(n_cells_to_plot)};
ax.XAxis.Visible = 'off';
ax.Color = 'none';
ax.YAxis.FontName ='arial';
ax.YAxis.FontSize = 8;
box off
ylabel('Neuron');
axis tight
set(gcf,'Units','centimeters');
set(gcf,'Position',[9 7 8 8]);

% Export
folder = 'C:\Users\mitch\OneDrive - University of Toronto\PhD\Presentations\Lab Meeting & Journal Club\LM6\Figures\Fig1';
print(gcf,fullfile(folder,'fig1e'),'-dpdf','-painters')
