function [sorted_rate_mtx_pref_dir, rate_mtx_pref_dir, pref_vector, rate_mtx_rightward, rate_mtx_leftward] = sort_cell_activity_prefDir(behavior_mtx, trace_mtx,trace_type, num_spatial_bins)
% plot rate matrix of spatial trial activity for every neuron in traces
% sort cells by location of FR peak
%
% trace mtx is (n,n,1)
%

% unique trials
[trials, ~, direction] = get_trials(behavior_mtx);

% travel direction
dir_vect = trial_direction(behavior_mtx);

% compute activity on every trial for every cell
%trace_type = 3;

rate_mtx_3d = binned_trl_activity_MDS(behavior_mtx, trace_mtx,trace_type, num_spatial_bins, trials);

% find preferred direction for each cell
rightward = direction == 0;
leftward = direction == 1;
active_right = sum((sum(rate_mtx_3d(rightward,:,:) > 0,2))>0); % Returns 1x1xneuron matrix with the number of trials which the cell is active
active_left = sum((sum(rate_mtx_3d(leftward,:,:) > 0,2))>0); % Repeat for left
fraction_active_right = squeeze(active_right)/sum(rightward); % Expressed as a fraction of trials
fraction_active_left = squeeze(active_left)/sum(leftward); % Repeat left
no_trial_activity = fraction_active_left == 0 & fraction_active_right == 0; % Identify cells with no trial activity
pref_right = fraction_active_right > fraction_active_left; % Logical index, marks cells more frequently active on rightward trials
pref_left = fraction_active_left > fraction_active_right; % Repeat for left
pref_right_cells = find(pref_right);
pref_left_cells = find(pref_left);

% Create a index vector where 1 = right prefering, 2 = left preferring, 3 =
% no trial activity
pref_vector = nan(size(pref_left));
pref_vector(pref_right) = 1;
pref_vector(pref_left) = 2;
pref_vector(no_trial_activity) = 3;


% calculate average place tuning for each direction
rate_mtx_rightward = squeeze(mean(rate_mtx_3d(rightward,:,:),1,'omitnan'))';
rate_mtx_leftward = squeeze(mean(rate_mtx_3d(leftward,:,:),1,'omitnan'))';


% average over preferred direction trials
rate_mtx_pref_dir = nan(size(rate_mtx_leftward));
rate_mtx_pref_dir(pref_right,:) = rate_mtx_rightward(pref_right,:);
rate_mtx_pref_dir(pref_left,:) = rate_mtx_leftward(pref_left,:);
sorted_rate_mtx_pref_dir = sort_rows_by_peak(rate_mtx_pref_dir);

% sort rows by peak
[sorted_rightward_right_pref, sort_idx_right_pref] = sort_rows_by_peak(rate_mtx_rightward(pref_right,:));
[sorted_leftward_left_pref, sort_idx_left_pref] = sort_rows_by_peak(rate_mtx_leftward(pref_left,:));


% plot sorted (normed)
figure;
tiledlayout(2,2,'TileSpacing','tight')
t1 = nexttile(1);
imagesc(sorted_rightward_right_pref);
title(sprintf('Right preferring cells \n rightward trials'));
set(gca,'YTick',1:size(sorted_rightward_right_pref,1),...
    'YTickLabel',cellstr(string(pref_right_cells(sort_idx_right_pref))));
ylabel('Cell Number')
xticklabels([])

t2 = nexttile(2);
imagesc(rate_mtx_leftward(pref_right_cells(sort_idx_right_pref),:));
title(sprintf('Right preferring cells \n leftward trials'));
set(gca,'YTickLabel',[])
ylabel([])
xticklabels([])

t3 = nexttile(3);
imagesc(sorted_leftward_left_pref);
title(sprintf('Left preferring cells \n leftward trials')) 
set(gca,'YTick',1:size(sorted_leftward_left_pref,1),...
    'YTickLabel',cellstr(string(pref_left_cells(sort_idx_left_pref))));
ylabel('Cell Number')
xticklabels([])

t4 = nexttile(4);
imagesc(rate_mtx_rightward(pref_left_cells(sort_idx_left_pref),:));
title(sprintf('Left preferring cells \n rightward trials'));
set(gca,'YTickLabel',[])
ylabel([])
xticklabels([])

clims = get([t1 t2 t3 t4],'clim');
set([t1 t2 t3 t4],'clim',max(cell2mat(clims)));





