
function generate_mouse_by_mouse_tuning_plots(group_all_cell_tcm_cell, group_pc_tcm_cell,group_names,mice_id_cell)

target_path = 'C:\Users\mitch\OneDrive - University of Toronto\PhD\Presentations\1on1_Meetings\2022-06-23\Mouse by Mouse Tuning Curve Correlations';

ngroups = length(group_all_cell_tcm_cell);

for igroup = 1:ngroups
    igroup_all_cell_tcm = group_all_cell_tcm_cell{igroup};
    igroup_pc_tcm = group_pc_tcm_cell{igroup};
    nmice = numel(igroup_all_cell_tcm);
    igroup_mice = mice_id_cell{igroup};
    for imouse = 1:nmice
        all_cell_tcm = igroup_all_cell_tcm{imouse};
        pc_tcm = igroup_pc_tcm{imouse};
        
        
        title_string_all = sprintf('%s mouse %s\nAll cells, no direction preference',group_names{igroup},igroup_mice{imouse});
        fname_all_cells = sprintf('%s_%s_all_cell',group_names{igroup},igroup_mice{imouse});
        days_apart_corr(all_cell_tcm,title_string_all);
        savefig(fullfile(target_path,fname_all_cells));
        print(gcf,fullfile(target_path, [fname_all_cells,'.svg']),'-dsvg','-painters');
        
        title_string_pc = sprintf('%s mouse %s\nPlace cells, preferred direction',group_names{igroup},igroup_mice{imouse});
        fname_pc_cells = sprintf('%s_%s_place_cell',group_names{igroup},igroup_mice{imouse});
        days_apart_corr(pc_tcm,title_string_pc);
        savefig(fullfile(target_path,fname_pc_cells));
        print(gcf,fullfile(target_path, [fname_pc_cells,'.svg']),'-dsvg','-painters');
    end
end


%{

        
        



for imouse = 1:numel(ctr_tcm_all)
    pc_mtx = ctr_pc_mtx{imouse};
    tcm_all = ctr_tcm_all{imouse};
    tcm_left = ctr_tcm_left{imouse};
    tcm_right = ctr_tcm_right{imouse};
    
    left_place_cells = any(pc_mtx == 1 | pc_mtx == 3,2);
    right_place_cells = any(pc_mtx == 2 | pc_mtx == 3,2);

    
    title_string_all = sprintf('%s mouse %s\nAll cells, no direction preference','control',ctr_mice{imouse});
    fp_all = sprintf('ctr mouse %s no dir pref',ctr_mice{imouse});
    days_apart_corr(tcm_all,title_string_all);
    savefig(fullfile(target_path,fp_all));
    print(gcf,fullfile(target_path, [fp_all,'.svg']),'-dsvg','-painters')
    
    title_string_left = sprintf('%s mouse %s\nLeft place cells, running left','control',ctr_mice{imouse});
    fp_left = sprintf('ctr mouse %s left',ctr_mice{imouse});
    days_apart_corr(tcm_left(left_place_cells,:,:),title_string_left);
    savefig(fullfile(target_path,fp_left));
    print(gcf,fullfile(target_path, [fp_left,'.svg']),'-dsvg','-painters')

    title_string_right = sprintf('%s mouse %s\nRight place cells, running_right','control',ctr_mice{imouse});
    fp_right = sprintf('ctr mouse %s right',ctr_mice{imouse});
    days_apart_corr(tcm_right(right_place_cells,:,:),title_string_right);
    savefig(fullfile(target_path,fp_right));
    print(gcf,fullfile(target_path, [fp_right,'.svg']),'-dsvg','-painters')
   
end

%%
for imouse = 1:numel(irr_tcm_all)
    pc_mtx = irr_pc_mtx{imouse};
    tcm_all = irr_tcm_all{imouse};
    tcm_left = irr_tcm_left{imouse};
    tcm_right = irr_tcm_right{imouse};
    
    left_place_cells = any(pc_mtx == 1 | pc_mtx == 3,2);
    right_place_cells = any(pc_mtx == 2 | pc_mtx == 3,2);

    
    title_string_all = sprintf('%s mouse %s\nAll cells, no direction preference','irradiation',irr_mice{imouse});
    fp_all = sprintf('irr mouse %s no dir pref',irr_mice{imouse});
    days_apart_corr(tcm_all,title_string_all);
    savefig(fullfile(target_path,fp_all));
    print(gcf,fullfile(target_path, [fp_all,'.svg']),'-dsvg','-painters')
    
    title_string_left = sprintf('%s mouse %s\nLeft place cells, running left','irradiation',irr_mice{imouse});
    fp_left = sprintf('irr mouse %s left',irr_mice{imouse});
    days_apart_corr(tcm_left(left_place_cells,:,:),title_string_left);
    savefig(fullfile(target_path,fp_left));
    print(gcf,fullfile(target_path, [fp_left,'.svg']),'-dsvg','-painters')

    title_string_right = sprintf('%s mouse %s\nRight place cells, running right','irradiation',irr_mice{imouse});
    fp_right = sprintf('irr mouse %s right',irr_mice{imouse});
    days_apart_corr(tcm_right(right_place_cells,:,:),title_string_right);
    savefig(fullfile(target_path,fp_right));
    print(gcf,fullfile(target_path, [fp_right,'.svg']),'-dsvg','-painters')
   
end

%%
for imouse = 1:numel(run_tcm_all)
    pc_mtx = run_pc_mtx{imouse};
    tcm_all = run_tcm_all{imouse};
    tcm_left = run_tcm_left{imouse};
    tcm_right = run_tcm_right{imouse};
    
    left_place_cells = any(pc_mtx == 1 | pc_mtx == 3,2);
    right_place_cells = any(pc_mtx == 2 | pc_mtx == 3,2);

    
    title_string_all = sprintf('%s mouse %s\nAll cells, no direction preference','exercise',run_mice{imouse});
    fp_all = sprintf('run mouse %s no dir pref',run_mice{imouse});
    days_apart_corr(tcm_all,title_string_all);
    savefig(fullfile(target_path,fp_all));
    print(gcf,fullfile(target_path, [fp_all,'.svg']),'-dsvg','-painters')
    
    title_string_left = sprintf('%s mouse %s\nLeft place cells, running left','exercise',run_mice{imouse});
    fp_left = sprintf('run mouse %s left',run_mice{imouse});
    days_apart_corr(tcm_left(left_place_cells,:,:),title_string_left);
    savefig(fullfile(target_path,fp_left));
    print(gcf,fullfile(target_path, [fp_left,'.svg']),'-dsvg','-painters')

    title_string_right = sprintf('%s mouse %s\nRight place cells, running right','exercise',run_mice{imouse});
    fp_right = sprintf('run mouse %s right',run_mice{imouse});
    days_apart_corr(tcm_right(right_place_cells,:,:),title_string_right);
    savefig(fullfile(target_path,fp_right));
    print(gcf,fullfile(target_path, [fp_right,'.svg']),'-dsvg','-painters')
   
end
%}