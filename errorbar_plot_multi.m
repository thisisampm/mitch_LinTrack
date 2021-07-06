function errorbar_plot_multi( cell_e , varargin)
%plots errorbar with different colors for the data points belonging to each
%row

% hold
hold on

% plot dots with colors
colors = distinguishable_colors(size(cell_e,1));
for igrp = 1:size(cell_e,1)
    errorbar_plot_noline(cell_e(igrp,:), [], [], colors(igrp,:))
end

% combine for overall errorbar
cell_e_comb = cell(1, size(cell_e,2));
for isamp = 1:size(cell_e,2)
    cell_e_comb{isamp} = [];
    for igrp = 1:size(cell_e,1)
        cell_e_comb{isamp} = [cell_e_comb{isamp}; cell_e{igrp,isamp}(:)];
    end
end

% plot overall errorbar
errorbar_plot_ebonly(cell_e_comb)