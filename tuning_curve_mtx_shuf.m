function [subj_means_obs, subj_means_shufs] = tuning_curve_mtx_shuf(tuning_curve_mtx, reference_session, num_shufs, varargin)
% computes the correlation between sessions for every neuron, as well as a
% correlation for every shuffle (# specificied by num_shufs) that randomizes 
% the rows on the first page
%
% tcm is a neurons X spaceBins X sessions matrix
%
% can input vector specifying which cells belong to which subject
%

%% inputs
if ~isempty(varargin)
    subj_grp = varargin{1};
else
    subj_grp = zeros(size(tuning_curve_mtx,1),1);
end
unq_subjs = unique(subj_grp);


%% compute observed correlation mean

% comparisons
comp_sessions = setdiff(1:size(tuning_curve_mtx,3), reference_session);

% subject means
subj_means_obs = nan(length(unq_subjs), length(comp_sessions));

% compute mean for each subject independently
for isubj = 1:length(unq_subjs)
    
    % local neurons
    subject_mtx = tuning_curve_mtx(subj_grp==unq_subjs(isubj),:,:);

    % local corrs
    subj_cell_corrs = subj_corrs(subject_mtx, comp_sessions, reference_session);
    
    % mean
    subj_means_obs(isubj,:) = nanmean(subj_cell_corrs);

end
    
    

%% shuffles

% subject means
subj_means_shufs = nan(length(unq_subjs), length(comp_sessions), num_shufs);

% iterate through shuffles
for ishuf = 1:num_shufs

    % compute mean for each subject independently
    for isubj = 1:length(unq_subjs)

        % local neurons
        subject_mtx = tuning_curve_mtx(subj_grp==unq_subjs(isubj),:,:);
        
        % shuffle reference session
        subject_mtx(:,:,reference_session) = subject_mtx(randperm(size(subject_mtx,1)),:,reference_session);

        % local corrs
        subj_cell_corrs = subj_corrs(subject_mtx, comp_sessions, reference_session);

        % mean
        subj_means_shufs(isubj,:,ishuf) = nanmean(subj_cell_corrs);

    end
end



%% plot
figure; hold on; 

% observed
ebp_in = cell(1,size(subj_means_obs,2));
for icomp = 1:length(ebp_in)
    ebp_in{icomp} = subj_means_obs(:,icomp);
end
for isubj = 1:length(unique(subj_grp))
    plot(subj_means_obs(isubj,:), 'o')
end
errorbar_plot_ebonly( ebp_in );

% shuffles (plot 95% center of overall means)
subj_means_shufs = squeeze(mean(subj_means_shufs,1))'; % shufs,comps
subj_means_shufs = sort(subj_means_shufs); % sorted
shuf_mean = mean(subj_means_shufs);
shuf_975 = subj_means_shufs(ceil(size(subj_means_shufs,1)*0.975),:);
shuf_025 = subj_means_shufs(ceil(size(subj_means_shufs,1)*0.025),:);
plot(shuf_mean, 'k-', 'linewidth', 2)
plot(shuf_975, 'k-', 'linewidth', 1)
plot(shuf_025, 'k-', 'linewidth', 1)

% aesthetics
ylim([-1 1])
hold on; plot(xlim, [0 0], 'k--')

end



%% internal functions
function subj_cell_corrs = subj_corrs(subject_mtx, comp_sessions, reference_session)

    % preallcate correlations
    subj_cell_corrs = nan(size(subject_mtx,1), size(subject_mtx,3)-1);

    % iterate through comparisons
    for icomp = 1:length(comp_sessions)

        % iterate through cells
        for ic = 1:size(subject_mtx,1)
            subj_cell_corrs(ic,icomp) = corr(subject_mtx(ic,:,reference_session)', subject_mtx(ic,:,comp_sessions(icomp))');
        end

    end

end




