%% Decoding figure


%% Perform decoding
% Across day decoding
[ctr_mean_error, ctr_median_error, ctr_true_positions, ctr_decoded_positions, ctr_mouse_ids] = decode_by_day1_for_group('reload\control');
[irr_mean_error, irr_median_error, irr_true_positions, irr_decoded_positions, irr_mouse_ids] = decode_by_day1_for_group('reload\irradiation');
[run_mean_error, run_median_error, run_true_positions, run_decoded_positions, run_mouse_ids] = decode_by_day1_for_group('reload\exercise');

%% Same day decoding
[ctr_sameday_mean_error, ctr_sameday_median_error, ctr_sameday_true_positions, ctr_sameday_decoded_positions, ctr_mouse_ids] = decode_by_sameday_for_group('reload\control');
[irr_sameday_mean_error, irr_sameday_median_error, irr_sameday_true_positions, irr_sameday_decoded_positions, irr_mouse_ids] = decode_by_sameday_for_group('reload\irradiation');
[run_sameday_mean_error, run_sameday_median_error, run_sameday_true_positions, run_sameday_decoded_positions, run_mouse_ids] = decode_by_sameday_for_group('reload\exercise');

%% Output folder
folder = 'C:\Users\mitch\OneDrive - University of Toronto\PhD\Presentations\Lab Meeting & Journal Club\LM6\Figures\Fig3';

%% Fig 3a. Population Vector example
% Get 4 tuning curves to stack to show population coding

load reload\control\196-3\2021-11-13.mat % Load an example session

rate_mtx_3d = binned_trl_activity_MDS(behavior_mtx,traces,3, 40, get_trials(behavior_mtx)); % Generate a trial x bin x neuron matrix
tuning_curves = squeeze(mean(rate_mtx_3d,1,'omitnan')); % Average tuning

example_cells = [34, 39, 114, 145];
t = tiledlayout(1,numel(example_cells),'TileSpacing','none','Padding','tight');

for i = 1:numel(example_cells)
    nexttile
    imagesc(tuning_curves(:,example_cells(i)));
    axis off
end

set(gcf,'Units','centimeters','Position',[5 5 3 12]);
%print(gcf,fullfile(folder,'tuning_curve_1d_examples'),'-dpdf','-painters')

%% Fig 3b. Decoding schematic

%% Fig 3c. Decoding error plots
color_mtx = [0 0 0; linspecer(2)];
group_names = {'ctr','irr','run'};
mouse_ids = {ctr_mouse_ids,irr_mouse_ids,run_mouse_ids};

% Plot cross-day decoding errors
day1_err = {ctr_median_error,irr_median_error,run_median_error};
sameday_err = {ctr_sameday_median_error,irr_sameday_median_error,run_sameday_median_error};

f = figure;
ax = axes;
hold(ax,'on');
linewidth = 1;
elapsed_time = 0:2:14;

for igroup = 1:numel(group_names)
    % Cross-day decoding
    igroup_day1_err = day1_err{igroup};
    mean_day1_err = mean(igroup_day1_err,2,'omitnan');
    std_day1_err = std(igroup_day1_err,[],2,'omitnan');
    n_mice = size(mouse_ids{igroup},1);
    sem_day1_err = std_day1_err./sqrt(n_mice);
    eb_day1(igroup) = errorbar(ax,elapsed_time(2:end),mean_day1_err,sem_day1_err,'LineWidth',linewidth,'Color',color_mtx(igroup,:),'CapSize',0);
    %legend(eb_day1,sprintf('%s Day 1 Decoder',group_names{igroup}))
end

for igroup = 1
 % Same-day decoding
    igroup_sameday_err = sameday_err{igroup};
    mean_sameday_err = mean(igroup_sameday_err,2,'omitnan');
    std_sameday_err = std(igroup_sameday_err,[],2,'omitnan');
    sem_sameday_err = std_sameday_err./sqrt(n_mice);
    eb_sameday(igroup) = errorbar(ax,elapsed_time,mean_sameday_err,sem_sameday_err,'LineWidth',linewidth,'Color',color_mtx(igroup,:),'CapSize',0,'LineStyle','--');
end
set(f,'Units','centimeters','Position',[5 5 10 10])
legend({'Control Day 1','Irradiation Day 1','Exercise Day 1','Sameday Decoding'},'FontName','arial','fontsize',6,'Location','best')
hold(ax,'off');
axis padded
ax.XAxis.FontName = 'arial';
ax.YAxis.FontName = 'arial';
ax.XAxis.FontSize = 8;
ax.YAxis.FontSize = 8;
ax.TickDir = 'out';
xlabel('Elapsed Time (days)','FontName','arial','FontSize',8);
ylabel('Median Decoder Error (cm)', 'FontName','arial','FontSize',8);
%legend(group_names)
%print(gcf,fullfile(folder,'Fig3_timelapse_decoding'),'-dpdf','-painters')




