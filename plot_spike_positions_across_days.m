function plot_spike_positions_across_days(mouse_folder,cell_reg_mtx,cell_to_plot, varargin)
% Loads the same cell trialwise activity on each day and plots

if ~isempty(varargin)
    running_only_indicator = varargin{1};
else
    running_only_indicator = 0;
end


files = get_file_paths_all(mouse_folder);
reg_id = contains(files,{'cell_reg','cellreg'},'IgnoreCase',1);
session_files = files(~reg_id);

figure;
t = tiledlayout(1,numel(session_files),'TileSpacing','tight','Padding','compact');

trace_type = 3;
nbins = 40;

cell_id = cell_reg_mtx(cell_to_plot,:);

for iday = 1:numel(session_files)
    load(session_files{iday},'behavior_mtx','traces');
    switch running_only_indicator
        case 1
            running_ts = behavior_mtx(:,5) > 5;
            behavior_mtx = behavior_mtx(running_ts,:);
            traces = traces(:,running_ts,:);
    end
        nexttile(iday)
    h = plot(behavior_mtx(:,2),behavior_mtx(:,1),'color',[0 0 0 0.5]); % Plot the mouse position
    axis padded
    if cell_id(iday) == 0
        title(sprintf('Inactive on Day %i',iday*2-1));
    else
        spike_ind = traces(cell_id(iday),:,1) == 1;
        hold on;
        plot(behavior_mtx(spike_ind,2), behavior_mtx(spike_ind,1),'r.','MarkerSize',15);
        hold off
        title(sprintf('Day %i',iday*2-1));
        axis padded
    end
end

sgtitle(sprintf('Cell %i',cell_to_plot));
    