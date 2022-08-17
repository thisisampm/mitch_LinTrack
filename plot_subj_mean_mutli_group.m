function [combined_group_table] = plot_subj_mean_mutli_group(group_tcm_cell,varargin)

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


%% For each group get subject mean tuning correlations

subj_mean_cell = cell(size(group_tcm_cell));



for igroup = 1:ngroups
    igroup_mice = mice_id_cell{igroup};
    nmice = numel(igroup_mice);
    group_tcm = group_tcm_cell{igroup}; % Pull out the group tuning curves from the larger cell
    nsessions = size(cat(1,group_tcm{:}),3);
    subj_means = nan(nsessions-1,nmice);
    group_index = repelem(group_names{igroup},nmice*(nsessions-1),1);
    mouse_index = [];
    session_separation = (1:(nsessions-1))';
    days_between_sessions = 2;
    days_apart_index = repmat(session_separation*days_between_sessions,nmice,1);
    
    for imouse = 1:nmice
        mouse_index = [mouse_index; repelem(igroup_mice(imouse),nsessions-1,1)];
        imouse_tcm = group_tcm{imouse};
        corr_mtx = days_apart_corr(imouse_tcm);      
        subj_means(:,imouse) = mean(corr_mtx,[1 3],'omitnan');
    end
    
    
    subj_mean_cell{igroup} = subj_means;
    group_table = table(group_index,mouse_index,days_apart_index, subj_means(:),...
            'VariableNames',{'group','mouse','days_apart','corr'});
    if igroup == 1
        combined_group_table = group_table;
    else
        combined_group_table = [combined_group_table;group_table];
    end
        
end

%% plot
color_mtx = [0 0 0; linspecer(ngroups-1)];

figure;hold on;
for igroup = 1:ngroups
    mouse_measurements = subj_mean_cell{igroup};
    mean_corr = mean(mouse_measurements,2,'omitnan');
    std_corr = std(mouse_measurements,[],2,'omitnan');
    n_sample = size(mouse_measurements,2);
    SEM = std_corr./sqrt(n_sample);
    upper_CI = mean_corr + tinv(0.975, n_sample -1).*SEM; % Calculate the upper 95% confidence interval for each
    lower_CI = mean_corr + tinv(0.025, n_sample -1).*SEM;
    day_separation = (0:2:14)';
    % Reverse Fisher Z transform for plotting
    mouse_measurements = [ones(1,size(mouse_measurements,2)); tanh(mouse_measurements)];
    mean_corr = [1;tanh(mean_corr)];
    upper_CI = [1;tanh(upper_CI)];
    lower_CI = [1;tanh(lower_CI)];
    SEM = [0; tanh(SEM)];
    for imouse = 1:n_sample
        plot(day_separation, mouse_measurements(:,imouse),'Color',color_mtx(igroup,:),'LineWidth',0.5);
    end
    plot(day_separation,(mean_corr),'LineWidth',4,'Color',color_mtx(igroup,:));
    patch([day_separation; flip(day_separation)],([upper_CI; flip(lower_CI)]),color_mtx(igroup,:),'edgecolor','none','facealpha',0.25);     
    %errorbar(day_separation,mean_corr,SEM,'LineWidth',4,'Color',color_mtx(igroup,:),'Capsize',0)
end
axis padded
legend(flip(findobj(gcf,'Linewidth',4)),group_names);
xlabel('Days Apart')
ylabel('Place Field Correlation (r)')











    
    
    






































