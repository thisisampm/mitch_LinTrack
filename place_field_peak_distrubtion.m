



function [place_field_peaks, fraction_pc] = place_field_peak_distrubtion(place_cell_info)

group_names = {'Control','Irradiation','Exercise'};
groups = unique(place_cell_info(:,1));
max_n_mice  = unique(place_cell_info(:,2));
nbins = numel(unique(place_cell_info(:,6)));

place_field_peaks = nan(numel(groups),numel(max_n_mice),nbins);
fraction_pc = nan(numel(groups),numel(max_n_mice));



for igroup = 1:numel(groups)
    group_index = place_cell_info(:,1) == igroup;
    mice = unique(place_cell_info(group_index,2));
    for imouse = 1:numel(mice)
        pc_id = place_cell_info(:,5) == 1;
        mouse_index = (place_cell_info(:,2) == imouse & group_index & pc_id);
        place_field_peaks(igroup,imouse,:) = histcounts(place_cell_info(mouse_index,6));
        m_index = (place_cell_info(:,2) == imouse);
        fraction_pc(igroup,imouse) = sum(m_index & group_index & pc_id)./sum(m_index & group_index);
    end
end


color_cell = mat2cell(linspecer(3),ones(1,3));
color_mtx = linspecer(3);

for i = 1:3
    mouse_id = squeeze(~isnan(place_field_peaks(i,:,:)));
    n_measurement(i) = sum(mouse_id(:,1));
end

mouse_pf_total = sum(place_field_peaks,3,'omitnan');
place_field_frac = place_field_peaks./mouse_pf_total;


pf_peak_means = squeeze(mean(place_field_frac,2,'omitnan'));
pf_peak_sem = squeeze(std(place_field_frac,[],2,'omitnan'))./sqrt(n_measurement');

figure;
eb = errorbar(pf_peak_means',pf_peak_sem');
set(eb,{'Color'}, color_cell,'LineWidth',1,'Capsize',0);
axis padded
title('Place Field Peak Distribution')
xlabel('Spatial Bin')
ylabel('Fraction of Place Field Peaks')
legend('Control','Irradiation','Exercise')

% fraction place cells
figure;
fraction_pc = fraction_pc*100;
fraction_pc_mean = mean(fraction_pc,2,'omitnan');
fraction_pc_sem = std(fraction_pc,[],2,'omitnan')./sqrt(n_measurement');
b1 = bar(fraction_pc_mean,'FaceColor','flat');
b1.CData = color_mtx;
xticklabels(group_names);
ylabel('Fraction of all cells that are place cells (%)')
hold on
err = errorbar(1:3,fraction_pc_mean, fraction_pc_sem);
err.Color = [0 0 0];
err.LineStyle = 'none';
err.CapSize = 0;
err.LineWidth = 1;
hold off


