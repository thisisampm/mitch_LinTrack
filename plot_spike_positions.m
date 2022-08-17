function [h] =  plot_spike_positions(behavior_mtx,traces, varargin)

% Takes behaviour matrix and calcium data from linear track imaging and
% plots the mouse position overlaid with a given cells spiking activity

if ~isempty(varargin)
    cell_id = varargin{1};
else
    cell_id = 1:size(traces,1);
end

ncells = numel(cell_id);
f = figure; % Open new figure for plotting

if ncells == 1
    ax = axes(f);
    h = plot(ax,behavior_mtx(:,2),behavior_mtx(:,1),'color',[0 0 0 0.2]); % Plot the mouse position
    spike_ind = traces(cell_id, :, 1) == 1; % Get an index of where the cell has spiked
    hold on;
    plot(behavior_mtx(spike_ind, 2), behavior_mtx(spike_ind, 1), 'r.','MarkerSize', 15); % Plot red dots at time and position of spikes
    hold off
    title(['Cell ' num2str(cell_id)]); % Title plot with cell number
    xlabel('Spatial Position');
    ylabel('Time');
elseif ncells >= 2 % If multiple cells subplot them
    t = tiledlayout(f,'flow');
    for i = 1:ncells
        nexttile(i)
        h = plot(behavior_mtx(:,2),behavior_mtx(:,1),'color',[0 0 0 0.2]); % Plot the mouse position
        spike_ind = traces(cell_id(i), :, 1) == 1; % Get an index of where the cell has spiked
        hold on;
        plot(behavior_mtx(spike_ind, 2), behavior_mtx(spike_ind, 1), 'r.','MarkerSize', 15); % Plot red dots at time and position of spikes
        hold off
        title(['Cell ' num2str(cell_id(i))]); % Title plot with cell number
        yticklabels([]);xticklabels([]);
    end
    t.Padding = 'tight';
    t.TileSpacing = 'compact';
    t.XLabel.String = 'Spatial Position';
    t.YLabel.String = 'Time';
end
        
    


