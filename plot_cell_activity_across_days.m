function plot_cell_activity_across_days(mouse_folder,cell_reg_mtx,cell_to_plot)
% Loads the same cell trialwise activity on each day and plots



files = get_file_paths_all(mouse_folder);
reg_id = contains(files,{'cell_reg','cellreg'},'IgnoreCase',1);
session_files = files(~reg_id);

figure;
tiledlayout(1,numel(session_files),'TileSpacing','tight','Padding','tight');

trace_type = 3;
nbins = 40;

cell_id = cell_reg_mtx(cell_to_plot,:);

for iday = 1:numel(session_files)
    load(session_files{iday},'behavior_mtx','traces');
    trials = get_trials(behavior_mtx);
    rate_mtx_3d = binned_trl_activity_MDS(behavior_mtx,traces,trace_type,nbins,trials);
    nexttile(iday);
    if cell_id(iday) == 0
        title(sprintf('Inactive on Session %i',iday));
    else
        imagesc(rate_mtx_3d(:,:,cell_id(iday)));
        title(sprintf('Day %i',iday*2-1));
    end
end

sgtitle(sprintf('Cell %i',cell_to_plot));
    
    