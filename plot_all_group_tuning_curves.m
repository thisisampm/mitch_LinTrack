function [] = plot_all_group_tuning_curves(all_tcm_irr, all_tcm_ctr, all_tcm_run)

ngroups = 3;
ndays = 8;
irr_pos = 1:8;
ctr_pos = 9:16;
run_pos = 17:24;
subplot_pos = [irr_pos; ctr_pos; run_pos];
%subplot_pos = 1:24; 

tcm_cell = {all_tcm_irr, all_tcm_ctr, all_tcm_run};

for igroup = 1:ngroups
    tcm = tcm_cell{igroup};
    for iday = 1:ndays
        subplot(ngroups,ndays,subplot_pos(igroup,iday));
        imagesc(norm_mtx(tcm(:,:,iday)')');
        yticklabels([])
        xticklabels([])
        if igroup == 1
           title_string = sprintf('Day %i',(iday*2-1));
           title(title_string);
        end
    end
end

%{
for i = 1:ndays
    subplot(1,ndays,i);
    imagesc(norm_mtx(tcm(:,:,i)')');
    yticklabels([])
    xticklabels([])
end
%}
    