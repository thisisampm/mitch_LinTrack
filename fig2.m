

% Figure 2
%%

[ctr_tcm_no_pref, ctr_tcm_left_pc, ctr_tcm_right_pc, ctr_pc_cell, ctr_pc_table, ctr_mice,ctr_pc_tcm] = extract_1D_spatial_tuning_multi('reload\control');
[irr_tcm_no_pref, irr_tcm_left_pc, irr_tcm_right_pc, irr_pc_cell, irr_pc_table, irr_mice,irr_pc_tcm] = extract_1D_spatial_tuning_multi('reload\irradiation');
[run_tcm_no_pref, run_tcm_left_pc, run_tcm_right_pc, run_pc_cell, run_pc_table, run_mice,run_pc_tcm] = extract_1D_spatial_tuning_multi('reload\exercise');
%%
folder = 'C:\Users\mitch\OneDrive - University of Toronto\PhD\Presentations\Lab Meeting & Journal Club\LM6\Figures\Fig2';
group_names = {'Control','Irradiation','Exercise'};
color_mtx = [0 0 0; linspecer(2)];

%% Fig 2a. Place coding example

%% Fig 2b. Pie chart showing the cells recorded from by group
%count_pc_stats % Run script count_pc_stats

active_means = mean(mean_p_active_by_session,2); % A group x session matrix with the fraction of cells active of total cells averaged by mouse

tiledlayout(2,3,'TileSpacing','tight','Padding','tight')
% Active/total cells
for i = 1:numel(group_names)
    nexttile
    pie([1- active_means(i) active_means(i)])
    title(group_names{i},'FontName','arial','FontSize',8)
    p = findobj(gca,'type','patch');
    txt = findobj(gca,'type','text');
    set(txt,'Fontname','arial','fontsize',8);
    set(p,'FaceColor',color_mtx(i,:));
    set(p(2),'FaceAlpha',0.25)
end

pc_both = squeeze(mean(n_pc_both./n_active_cells,[2 3],'omitnan'));
pc_left = squeeze(mean(n_pc_left./n_active_cells,[2 3],'omitnan'));
pc_right = squeeze(mean(n_pc_right./n_active_cells,[2 3],'omitnan'));
pc_none = 1 - pc_left - pc_right - pc_both;

%Place coding/active
pc_color_mtx = linspecer(6);
pc_color_mtx = pc_color_mtx(3:end,:);
pc_color_mtx = distinguishable_colors(4);
for i = 1:numel(group_names)
    nexttile
    pie([pc_none(i) pc_left(i) pc_both(i) pc_right(i)])
    %title(group_names{i},'FontName','arial','FontSize',8)
    p = findobj(gca,'type','patch');
    txt = findobj(gca,'type','text');
    set(txt,'Fontname','arial','fontsize',8);
    set(p,{'facecolor'},num2cell(pc_color_mtx,2))
    set(p,'FaceAlpha',0.65)
end

set(gcf,'Units','centimeters','Position',[5 5 8 5]);

%print(gcf,fullfile(folder,'place_cell_proportions'),'-dpdf','-painters')


%%


% determine cells that are not active on day 1 and only plot active place
% cells

group_pc_tcm = {irr_pc_tcm,ctr_pc_tcm,run_pc_tcm};
group_names = {'Irradiation','Control','Exercise'};
n_sessions = size(ctr_pc_tcm{1},3);
n_groups = numel(group_names);

color_mtx = color_mtx([2 1 3],:); % Rearrange for this plot

t = tiledlayout(n_groups,n_sessions,'TileSpacing','tight','Padding','tight');

tile_count = 1;
for igroup = 1:numel(group_names)
    igroup_tcm = group_pc_tcm{igroup};
    
    combined_pc_tcm = cat(1,igroup_tcm{:});
    day1_inactive = (all(isnan(combined_pc_tcm(:,:,1)),2) | all(combined_pc_tcm(:,:,1)==0,2));
    combined_pc_tcm = combined_pc_tcm(~day1_inactive,:,:); % Restrict to active cells for plotting
    [~, day1_sorting] = sort_rows_by_peak(combined_pc_tcm(:,:,1)); % Index sorting day1 field peaks

    n_sessions = size(combined_pc_tcm,3);
    for isesh = (1:n_sessions)
        nexttile(tile_count)
        tile_count = tile_count + 1;
        imagesc((combined_pc_tcm(day1_sorting,:,isesh)')');
        %colormap gray
        caxis([0 1]);
        % Session labels on top
        if igroup == 1
            title(sprintf('Day %i',isesh*2-1),'FontName','arial','FontSize',8);
        end
        % Group labels on left
        if isesh == 1
            ylabel(group_names{igroup},'FontName','Arial','FontSize',8,'Color',color_mtx(igroup,:));
            set(get(gca,'yaxis'),'fontname','arial','fontsize',8)
            xticklabels([]);
            box off
            set(gca,'TickLength',[0 0])
        else
            axis off
        end
        % One colorbar for whole group
        if igroup == 2 && isesh == 8
            cb = colorbar('FontName','arial','FontSize',8,'Ticks',[0 1]);
            cb.Label.String = 'p(Active|Position)';
        end
    end
end
color_mtx = [0 0 0; linspecer(2)];

set(gcf,'Units','centimeters','Position',[5 5 15 10]);
%print(gcf,fullfile(folder,'all_group_all_session_tuning_curves'),'-dpdf','-painters')

%% Fig 2d. Activity map correlation example
t2 = tiledlayout(2,2,'TileSpacing','compact','Padding','tight');
nexttile
%example_cell = 15;
example_cell = 1;
% E.g. high correlation
plot(smooth(combined_pc_tcm(example_cell,:,1),5),'LineWidth',2,'color',[0 0 0 0.5])
ylabel('Activity','FontName','arial','FontSize',8)
title('Day i','FontName','arial','FontSize',8)
axis padded
ylim_hold = get(gca,'ylim');
set(gca,'YTick',[],'XTick',[]);

nexttile
plot(smooth(combined_pc_tcm(example_cell,:,2),5),'LineWidth',2,'color',[0 0 0 0.5])
title('Day i + n','FontName','arial','FontSize',8)
axis padded
set(gca,'ylim',ylim_hold)
set(gca,'YTick',[],'XTick',[]);
% E.g. low correlation
example_cell_lowr = 322;
nexttile
plot(smooth(combined_pc_tcm(example_cell_lowr,:,1),5),'LineWidth',2,'color',[0 0 0 0.5])
ylabel('Activity','FontName','arial','FontSize',8)
%title('Day i','FontName','arial','FontSize',8)
axis padded
ylim_hold = get(gca,'ylim');
set(gca,'YTick',[],'XTick',[]);

nexttile
plot(smooth(combined_pc_tcm(example_cell_lowr,:,5),5),'LineWidth',2,'color',[0 0 0 0.5])
%title('Day i + n','FontName','arial','FontSize',8)
axis padded
set(gca,'ylim',ylim_hold)
set(gca,'YTick',[],'XTick',[]);


sgtitle('Spatial Position','Fontname','arial','fontsize',8)
set(gcf,'Units','centimeters','Position',[5 5 6 5])
%print(gcf,fullfile(folder,'tuning_curve_example'),'-dpdf','-painters')



%% Fig 2e. Days apart correlation calculation
[combined_all_cell_corr_table] = days_apart_corr_multi({ctr_pc_tcm,irr_pc_tcm,run_pc_tcm},{'ctr','irr','run'},{ctr_mice,irr_mice,run_mice});
f2 = gcf;
set(f2.Children,'Fontname','arial','fontsize',8)
set(f2,'Units','centimeters','Position',[5 5 6 10])
set(gca,'TickDir','out')
%print(gcf,fullfile(folder,'all_cell_pf_corr'),'-dpdf','-painters')

% Stats
% Fit a linear mixed effects model accounting for the effect of 
corr_lme = fitlme(combined_all_cell_corr_table,'corr ~ group + days_apart + group*days_apart +(1|mouse)+(1|mouse:cell)')





























