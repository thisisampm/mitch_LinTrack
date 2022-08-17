function [combined_corr_table, corr_mtx_cell] = days_apart_corr_multi(group_tcm_cell,varargin)

% group_tcm_cell should be a cell that contains a cell array with the
% tuning matrices for each mouse in that group

ngroups = length(group_tcm_cell); % The number of groups corresponds to the number of cells in the outermost cell array

% Optionally take the group names as a cell array the same size as group
% tcm.
% Optionally take the names of the mice in each group as cell's in an array
% the size of group tcm
if length(varargin) == 1
    group_names = varargin{1};
     nmice_per_group = cellfun(@length, group_tcm_cell);
    mice_id_cell = cell(size(nmice_per_group));
    for igroup = 1:length(mice_id_cell)
        mice_id_cell{igroup} = 1:nmice_per_group(igroup);
    end
elseif length(varargin) == 2
    group_names = varargin{1};
    mice_id_cell = varargin{2};
else
    group_names = num2cell(1:ngroups); % If no group names given, number them 1 through n groups
    nmice_per_group = cellfun(@length, group_tcm_cell);
    mice_id_cell = cell(size(nmice_per_group));
    for igroup = 1:length(mice_id_cell)
        mice_id_cell{igroup} = 1:nmice_per_group(igroup);
    end
end


corr_mtx_cell = cell(size(group_tcm_cell));
subj_index_cell = cell(size(group_tcm_cell));
group_index_cell = cell(size(group_tcm_cell));


% For each group combine the cells from each subject and create an index
% that identifies which mouse contributed which tuning curves
cell_count = 0;
for igroup = 1:ngroups
    subj_cell = group_tcm_cell{igroup}; % Open the current group's cell of tuning curves
    nsubj = length(subj_cell);
    ncell_per_subj = cellfun('size',subj_cell,1); % For each subject, determine the number of cells as size of first array dimension
    % Construct subject index vector
    subj_index = [];
    igroup_mice = mice_id_cell{igroup};
    cell_index = [];
    for isubj = 1:nsubj
        subj_index = [subj_index; repelem(igroup_mice(isubj),ncell_per_subj(isubj),1)];
        %cell_index = [cell_index;(1:ncell_per_subj(isubj))'];
        cell_index = [cell_index; (cell_count+1:1:ncell_per_subj(isubj)+cell_count)'];
        cell_count = cell_index(end);
    end
    subj_index_cell{igroup} = subj_index;
    all_group_tcm = cat(1,subj_cell{:}); % Concatenate all tuning curves from all mice in a group
    group_corr_mtx = days_apart_corr(all_group_tcm); % Calculate correlation of same cell's tuning to itself on different days
    group_index = repelem(group_names{igroup},length(group_corr_mtx),1);
    group_index_cell{igroup} = group_index;
    corr_mtx_cell{igroup} = group_corr_mtx;
    
    ncells = size(group_corr_mtx,1);
    n_session_separation = size(group_corr_mtx,2);
    n_session_compared = size(group_corr_mtx,3);
    
    % Take the 3d dimension of corr_mtx, corresponding to which sessions
    % are used and stack them under the columns (corresponding to how far
    % apart the sessions used in the correlation are)
    group_corr_mtx_2d = num2cell(group_corr_mtx,[1 2]);
    group_corr_mtx_2d = cat(1,group_corr_mtx_2d{:});
    group_corr_mtx_1d = group_corr_mtx_2d(:);
    % Index that indicates how far apart the two sessions used in the corr
    days_between_session = 2;
    day_separation = (1:n_session_separation)*days_between_session;
    session_separation_index = repmat(day_separation,n_session_separation*ncells,1);
    session_separation_index = session_separation_index(:); % Linearize
    
    % Indices are currently the length of dim1 of 3d_corr matrix, need to
    % replicate to match linearized length
    group_index = (cellstr(repmat(group_index,n_session_compared*n_session_separation,1)));
    subj_index = (cellstr(repmat(subj_index,n_session_compared*n_session_separation,1)));
    cell_index = repmat(cell_index,n_session_compared*n_session_separation,1);
    
    group_table = table(group_index, subj_index, cell_index, session_separation_index, group_corr_mtx_1d,...
        'VariableNames',{'group','mouse','cell','days_apart','corr'});
    if igroup == 1
        combined_corr_table = group_table;
    else
        combined_corr_table = [combined_corr_table; group_table];
    end
end

% Remove all the rows that do not have a correlation measurement
combined_corr_table = combined_corr_table(~isnan(combined_corr_table.corr),:);

%% Plot
color_mtx = [0 0 0; linspecer(ngroups-1)];

figure; hold on;

for igroup = 1:ngroups
    % Calculate summary statistics for each group
    mean_days_apart_corr = mean(corr_mtx_cell{igroup},[1 3],'omitnan'); % Mean Fisher Z correlation for all cells n days apart
    std_days_apart_corr = std(corr_mtx_cell{igroup},[],[1 3],'omitnan'); % Standard deviation
    n_sample_days_aprt_corr = sum(~isnan(corr_mtx_cell{igroup}),[1 3]); % The number of correlation measurements at each day separation
    SEM = std_days_apart_corr./sqrt(n_sample_days_aprt_corr); % The standard error of the mean
    upper_CI = mean_days_apart_corr + tinv(0.975, n_sample_days_aprt_corr -1).*SEM; % Calculate the upper 95% confidence interval for each
    lower_CI = mean_days_apart_corr + tinv(0.025, n_sample_days_aprt_corr -1).*SEM;
    
    % Plot swarmchart of correlations on each day
    for icomp = 1:numel(day_separation)
        icomp_corr = corr_mtx_cell{igroup}(:,icomp,:);
        swarmchart(day_separation(icomp)*ones(numel(icomp_corr),1),icomp_corr(:),...
            'markeredgealpha',0.01,'markeredgecolor',color_mtx(igroup,:));
    end
    % Plot mean plus 95% confidence interval
    plot(day_separation,mean_days_apart_corr,'LineWidth',3,'Color',color_mtx(igroup,:));
    patch([day_separation fliplr(day_separation)],[upper_CI fliplr(lower_CI)],color_mtx(igroup,:),'edgecolor','none','facealpha',0.5);     
end

hold off;
    maxr = .99; 
    tic_vect = [-maxr -.9 -.8 -.6 -.3 0 .3 .6 .8 .9 maxr];
    yticks(atanh(tic_vect));
    yticklabels(tic_vect);
    xticks(day_separation);
    low_r = min(cat(1,corr_mtx_cell{:}),[],'all');
    ylim([low_r,atanh(maxr)]);
    xlabel('Days Apart')
    ylabel('Place Field Correlation (r)')
    legend(flip(findobj(gcf,'type','line')),(group_names))


    













    
        
    
