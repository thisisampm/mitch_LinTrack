%function []

group_pc_cells = {ctr_pc_cell; irr_pc_cell; run_pc_cell};
ngroup = numel(group_pc_cells);
n_mice_per_group = cellfun('size',group_pc_cells,1); % Each mouse has its own cell so get size of cell array for n mice
group_names = {'ctr','irr','run'};

total_unique_cells = nan(size(group_pc_cells)); % Count of all unique cells in a group

n_sessions = cellfun('size',group_pc_cells{1},2);
n_sessions = n_sessions(1);
n_active_cells = nan(ngroup, max(n_mice_per_group),n_sessions);
n_unique_cells_per_mouse = nan(ngroup, max(n_mice_per_group),n_sessions);
n_pc_left = nan(ngroup,max(n_mice_per_group),n_sessions);
n_pc_right = nan(ngroup,max(n_mice_per_group),n_sessions);
n_pc_both = nan(ngroup,max(n_mice_per_group),n_sessions);
n_pc_any = nan(ngroup,max(n_mice_per_group),n_sessions);



for igroup = 1:ngroup
    igroup_pc_cell = group_pc_cells{igroup};
    total_unique_cells(igroup) = size(cat(1,igroup_pc_cell{:}),1); % combine all rows representing unique cells and count
    
    nmice = numel(igroup_pc_cell);
    %n_unique_cells_per_mouse = nan(1,nmice);
    %n_sessions = size(igroup_pc_cell{1},2);
    %n_active_cells_by_day = nan(nmice,n_sessions);
    
    for imouse = 1:nmice
        mouse_pc_mtx = igroup_pc_cell{imouse};
        n_unique_cells_per_mouse(igroup,imouse,:) = size(mouse_pc_mtx,1);
        n_active_cells(igroup,imouse,:) = sum(~isnan(mouse_pc_mtx),1);
        %n_active_cells_by_day(imouse, :) = sum(~isnan(mouse_pc_mtx),1);
        n_pc_left(igroup,imouse,:) = sum(mouse_pc_mtx == 1,1);
        n_pc_right(igroup,imouse,:)= sum(mouse_pc_mtx == 2,1);
        n_pc_both(igroup,imouse,:) = sum(mouse_pc_mtx == 3,1);
        n_pc_any(igroup,imouse,:) = sum(mouse_pc_mtx == 1 | mouse_pc_mtx == 2 | mouse_pc_mtx == 3,1);
    end
end
%% Summary stats
p_active_cells = n_active_cells./n_unique_cells_per_mouse;
p_pc_left = n_pc_left./n_active_cells;
p_pc_right= n_pc_right./n_active_cells;
p_pc_both = n_pc_both./n_active_cells;
p_pc_any = n_pc_any./n_active_cells;

%fraction active
mean_p_active_by_session = squeeze(mean(p_active_cells,2,'omitnan'));
std_p_active_by_session = squeeze(std(p_active_cells,[],2,'omitnan'));
SEM_active = std_p_active_by_session./sqrt(n_mice_per_group);
%fraction place coding
mean_pc_any_by_session = squeeze(mean(p_pc_any,2,'omitnan'));
std_pc_any_by_session = squeeze(std(p_pc_any,[],2,'omitnan'));
SEM_pc_any = std_p_active_by_session./sqrt(n_mice_per_group);

%% plot
%fraction active
color_mtx = [0 0 0; linspecer(2)];
figure;tiledlayout(2,1);nexttile;hold on
for igroup = 1:size(SEM_active,1)
    errorbar(mean_p_active_by_session(igroup,:),SEM_active(igroup,:),'LineWidth',1,'Color',color_mtx(igroup,:),'CapSize',0)
end
hold off
legend(group_names)
axis padded
ylim([0 1])
title('Fraction of total cells active')
ylabel('p(active|total)')
xticklabels(cellstr(string(0:2:14)));
xlabel('Day')
set(gca,'TickDir','out')

% fraction place coding
nexttile;  hold on
for igroup = 1:size(SEM_active,1)
    errorbar(mean_pc_any_by_session(igroup,:),SEM_pc_any(igroup,:),'LineWidth',1,'Color',color_mtx(igroup,:),'CapSize',0)
end
hold off
legend(group_names)
axis padded
ylim([0 1])
title('Fraction of active cells that are place cells')
ylabel('p(place coding|active)')
xticklabels(cellstr(string(0:2:14)));
xlabel('Day')
set(gca,'TickDir','out')

% Update for 
set(findall(gcf,'-property','FontSize'),'FontSize',8,'FontName','arial')
set(gcf,'Units','centimeters','Position',[5 5 6.5 10])




%% Place field peaks
ctr_pf_peak_width = [ctr_pc_table.width_left(ctr_pc_table.Place_cell ==1|ctr_pc_table.Place_cell ==3);...
    ctr_pc_table.width_right(ctr_pc_table.Place_cell == 2 | ctr_pc_table.Place_cell == 3)];

irr_pf_peak_width = [irr_pc_table.width_left(irr_pc_table.Place_cell ==1|irr_pc_table.Place_cell ==3);...
    irr_pc_table.width_right(irr_pc_table.Place_cell == 2 | irr_pc_table.Place_cell == 3)];

run_pf_peak_width = [run_pc_table.width_left(run_pc_table.Place_cell ==1|run_pc_table.Place_cell ==3);...
    run_pc_table.width_right(run_pc_table.Place_cell == 2 | run_pc_table.Place_cell == 3)];

ecdf(ctr_pf_peak_width)
hold on
ecdf(irr_pf_peak_width)
ecdf(run_pf_peak_width)
hold off

%%
mice_cell = {ctr_mice,irr_mice,run_mice};
days = 0:2:14;
output_folder = 'C:\Users\mitch\OneDrive - University of Toronto\PhD\Presentations\1on1_Meetings\2022-06-23\Mouse by Mouse Active Fraction';
for igroup = 1:numel(mice_cell)
    mice = mice_cell{igroup};
    for imouse = 1:numel(mice)
        figure;
        plot(days, squeeze(n_active_cells(igroup,imouse,:)),'LineWidth',3);
        hold on
        plot(days, squeeze(n_pc_any(igroup,imouse,:)),'LineWidth',3);
        hold off
        title_string = sprintf('%s mouse %s',group_names{igroup},mice{imouse});
        title(title_string);
        axis padded
        ylim([0 270])
        legend('Active fraction','Place cell fraction')
        ylabel('n cells')
        xlabel('Day')
        %savefig(fullfile(output_folder,title_string));
        %print(gcf,fullfile(output_folder, [title_string,'.svg']),'-dsvg','-painters');
    end
end
    
    





















