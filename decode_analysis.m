%% Decoding analysis driver script for lab meeting 5
% 2022-04-19

if ~exist('ctr_median_error','var')
    load('2022-04-19_Decoding_workspace.mat')
end

color_mtx = linspecer(3); % 3 distinguishable colors to use for each group
group_id = {'Control','Irradiation','Exercise'};

%% Errorbar comparing

cross_day_x = 2:2:14; % Number of days decoder is away from training data for cross day decoding
sameday_x = 0:2:14; % vectors of days with same day recording

% Control timelapse vs same day decoding
subplot_f = figure;
ctr_crossday_means = mean(ctr_median_error,2);
ctr_crossday_sem = std(ctr_median_error,[],2)./size(ctr_median_error,2);
ctr_ax = subplot(1,3,1);
ctr_h1 = errorbar(ctr_ax, cross_day_x, ctr_crossday_means,ctr_crossday_sem);
set(ctr_h1, 'LineWidth', 2, 'Capsize', 0, 'Color', color_mtx(1,:))
hold on
ctr_sameday_means = mean(ctr_sameday_median_error,2);
ctr_sameday_sem = std(ctr_sameday_median_error,[],2)./size(ctr_median_error,2);
ctr_h2 = errorbar(ctr_ax, sameday_x,ctr_sameday_means,ctr_sameday_sem);
hold off
set(ctr_h2, 'LineWidth', 2, 'Capsize', 0, 'Color', [0.5 0.5 0.5])
axis padded
xlabel('Elapsed Time (Days)')
ylabel('Median Decoder Error (Spatial bins)')
title('Control')
legend('Day 1 Decoder','Same Day Decoder','Location','northwest')

% Irradiation timelapse vs sameday decoding
irr_crossday_means = mean(irr_median_error,2);
irr_crossday_sem = std(irr_median_error,[],2)./size(irr_median_error,2);
irr_ax = subplot(1,3,2);
irr_h1 = errorbar(irr_ax, cross_day_x, irr_crossday_means,irr_crossday_sem);
set(irr_h1, 'LineWidth', 2, 'Capsize', 0, 'Color', color_mtx(2,:))
hold on
irr_sameday_means = mean(irr_sameday_median_error,2);
irr_sameday_sem = std(irr_sameday_median_error,[],2)./size(irr_median_error,2);
irr_h2 = errorbar(irr_ax,sameday_x,irr_sameday_means,irr_sameday_sem);
hold off
set(irr_h2, 'LineWidth', 2, 'Capsize', 0, 'Color', [0.5 0.5 0.5])
axis padded
xlabel('Elapsed Time (Days)')
ylabel('Median Decoder Error (Spatial bins)')
title('Irradiation')
legend('Day 1 Decoder','Same Day Decoder','Location','northwest')

% Exercise Timelapse vs same day decoding
run_crossday_means = mean(run_median_error,2);
run_crossday_sem = std(run_median_error,[],2)./size(run_median_error,2);
run_ax = subplot(1,3,3);
run_h1 = errorbar(run_ax, cross_day_x, run_crossday_means,run_crossday_sem);
set(run_h1, 'LineWidth', 2, 'Capsize', 0, 'Color', color_mtx(3,:))
hold on
run_sameday_means = mean(run_sameday_median_error,2);
run_sameday_sem = std(run_sameday_median_error,[],2)./size(run_median_error,2);
run_h2 = errorbar(run_ax,sameday_x,run_sameday_means,run_sameday_sem);
hold off
set(run_h2, 'LineWidth', 2, 'Capsize', 0, 'Color', [0.5 0.5 0.5])
axis padded
xlabel('Elapsed Time (Days)')
ylabel('Median Decoder Error (Spatial bins)')
title('Exercise')
legend('Day 1 Decoder','Same Day Decoder','Location','northwest')

ylims = get([ctr_ax, irr_ax, run_ax],'ylim');
ylims = cell2mat(ylims);
ylim_upper = max(ylims(:));
set([ctr_ax irr_ax run_ax],'ylim',[0 ylim_upper]);

% 3 Group Timelapse
all_f = figure;
all_ax = axes;
ctr_h3 = errorbar(all_ax, cross_day_x, ctr_crossday_means,ctr_crossday_sem);
set(ctr_h3, 'LineWidth', 2, 'Capsize', 0, 'Color', color_mtx(1,:))
hold on
irr_h3 = errorbar(all_ax, cross_day_x, irr_crossday_means,irr_crossday_sem);
set(irr_h3, 'LineWidth', 2, 'Capsize', 0, 'Color', color_mtx(2,:))
run_h3 = errorbar(all_ax, cross_day_x, run_crossday_means,run_crossday_sem);
set(run_h3, 'LineWidth', 2, 'Capsize', 0, 'Color', color_mtx(3,:))
hold off
xlabel('Elapsed Time (Days)')
ylabel('Median Decoder Error (Spatial bins)')
axis padded
legend(group_id, 'Location','northwest')

%% Repeat with dropping sessions with poor sameday decoding (err > 6)
error_threshold = 6;
ctr_high_error = ctr_sameday_median_error > error_threshold;
irr_high_error = irr_sameday_median_error > error_threshold;
run_high_error = run_sameday_median_error > error_threshold;
ctr_high_error = ctr_high_error(2:end,:); % Cut off top row corresponding day 1
irr_high_error = irr_high_error(2:end,:); % Cut off top row corresponding day 1
run_high_error = run_high_error(2:end,:); % Cut off top row corresponding day 1

ctr_median_error2 = ctr_median_error;
ctr_median_error2(ctr_high_error) = NaN;


irr_median_error2 = irr_median_error;
irr_median_error2(irr_high_error) = NaN;


run_median_error2 = run_median_error;
run_median_error2(run_high_error) = NaN;

ctr_crossday_means2 = mean(ctr_median_error2,2,'omitnan');
ctr_crossday_sem2 = std(ctr_median_error2,[],2,'omitnan')./size(ctr_median_error2,2);

irr_crossday_means2 = mean(irr_median_error2,2,'omitnan');
irr_crossday_sem2 = std(irr_median_error2,[],2,'omitnan')./size(irr_median_error2,2);

run_crossday_means2 = mean(run_median_error2,2,'omitnan');
run_crossday_sem2 = std(run_median_error2,[],2,'omitnan')./size(run_median_error2,2);

% 3 Group Timelapse
all_f2 = figure;
all_ax2 = axes;
ctr_h4 = errorbar(all_ax2, cross_day_x, ctr_crossday_means2,ctr_crossday_sem2);
set(ctr_h4, 'LineWidth', 2, 'Capsize', 0, 'Color', color_mtx(1,:))
hold on
irr_h4 = errorbar(all_ax2, cross_day_x, irr_crossday_means2,irr_crossday_sem2);
set(irr_h4, 'LineWidth', 2, 'Capsize', 0, 'Color', color_mtx(2,:))
run_h4 = errorbar(all_ax2, cross_day_x, run_crossday_means2,run_crossday_sem2);
set(run_h4, 'LineWidth', 2, 'Capsize', 0, 'Color', color_mtx(3,:))
hold off
xlabel('Elapsed Time (Days)')
ylabel('Median Decoder Error (Spatial bins)')
axis padded
legend(group_id, 'Location','northwest')


