
target_path = 'C:\Users\mitch\OneDrive - University of Toronto\PhD\Presentations\1on1_Meetings\2022-06-23\Mouse by Mouse Tuning Curve Correlations';




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