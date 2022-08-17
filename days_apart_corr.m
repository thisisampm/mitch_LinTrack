%% Function that compares the tuning curves of all cells n days apart

% Inputs
% all_tcm: a neuron x spatial bin x session matrix with neurons ordered
% according to Cell Reg cell_to_index_map. (2nd output from
% place_cell_tuning.)
% Outputs
% A cell's correlation to itself

function [corr_mtx, session_compared_index] = days_apart_corr(all_tcm, varargin)


if length(varargin) == 1
    plot_cond = true;
    title_string = string(varargin{1});
elseif length(varargin) == 2
    plot_cond = true;
    title_string = string(varargin{1});
    tile_target = varargin{2};
else 
    plot_cond = false;
end
    

n_cells = size(all_tcm,1);
n_bins = size(all_tcm,2);
n_sessions = size(all_tcm,3);

% Identify all possible n+1,n+3,n+5...n+max_separation comparisons that can be made
sessions = 1:n_sessions; % Vector of all sessions
n_sessions_apart_mtx = sessions - sessions'; % Matrix indicating how far each session is from every other session

session_separations = unique(n_sessions_apart_mtx(n_sessions_apart_mtx > 0)); % Vector of the unique sesssion separation possibilities i.e. n+1=1, n+2=2, ... n+max sepeartion
n_separations = numel(session_separations);
n_comparisons = sum(n_sessions_apart_mtx > 0,[1 2]); % Sum the logical vector along both dimensions to get the number of all comparisons
[~, max_n_comparisons] = mode(n_sessions_apart_mtx(n_sessions_apart_mtx>0),[1 2]); % The number of times which the most frequent comparison distance appears

corr_mtx = nan(n_cells,n_separations, max_n_comparisons); % Preallocate a cell x seperation distance x comparison output matrix

% For each possible comparison calculate the correlation between a cell's
% tuning curve
comp_count = 0;
session_compared_index = nan(n_comparisons, 2); % Column vector indicating which two sessions are compared for each column of output corr matrix
for i_separation = session_separations'
    [session1, session2] = find(n_sessions_apart_mtx == i_separation); % From the sessions apart matrix, determine which sessions can be compared at each separation
    for icomp = 1:length(session1)
        comp_count = comp_count + 1;
        sesh_to_sesh_corr = diag(corr(all_tcm(:,:,session1(icomp))',all_tcm(:,:,session2(icomp))')); % Compare tuning curve of each cell on a given session to itself in the other
        sesh_to_sesh_corr(sesh_to_sesh_corr > 0.9999) = 0.9999; % Restrict correlations to 6 significant digits because correlations too close to 1 will produce Inf in Fisher Z transform
        corr_mtx(:,i_separation,icomp) = atanh(sesh_to_sesh_corr); % Fisher Z transform the Pearson correlation output
        
        % Keep track of where each comparison comes from
        session_compared_index(comp_count,1) = session1(icomp);
        session_compared_index(comp_count,2) = session2(icomp);
    end
end

%% Calculate a chance distribution by shuffling cell identities

nshuffles = 500;
for ishuf = 1:nshuffles
end
    

%% Calculate summary statistics

mean_days_apart_corr = mean(corr_mtx,[1 3],'omitnan'); % Mean Fisher Z correlation for all cells n days apart
std_days_apart_corr = std(corr_mtx,[],[1 3],'omitnan'); % Standard deviation
n_sample_days_aprt_corr = sum(~isnan(corr_mtx),[1 3]); % The number of correlation measurements at each day separation
SEM = std_days_apart_corr./sqrt(n_sample_days_aprt_corr); % The standard error of the mean
upper_CI = mean_days_apart_corr + tinv(0.975, n_sample_days_aprt_corr -1).*SEM; % Calculate the upper 95% confidence interval for each
lower_CI = mean_days_apart_corr + tinv(0.025, n_sample_days_aprt_corr -1).*SEM;


%% Plot
if plot_cond

    days_between_sessions = 2;
    day_separation = session_separations*days_between_sessions;
    
    % If passed a tiledlayout object, plot to that. Otherwise new fig
    if exist('tile_target','var')
        nexttile(tile_target);
    else
        figure;
    end
    
    hold on;
    for icomp = 1:numel(day_separation)
        icomp_corr = corr_mtx(:,icomp,:);
        swarmchart(day_separation(icomp)*ones(numel(icomp_corr),1),icomp_corr(:),'MarkerEdgeAlpha',0.2,'MarkerEdgeColor',[0 0 0]);
    end
    l = plot(day_separation, mean_days_apart_corr,'LineWidth',3,'Color',[0 0 0]);
    p = patch([day_separation; flip(day_separation)],[upper_CI'; flip(lower_CI')],[1 0 0],'EdgeColor','none','FaceAlpha',0.25);
    legend([l,p],'Mean','95% CI');
    %plot(day_separation, upper_CI,'LineWidth',2,'LineStyle','--','Color',[0 0 0]);
    %plot(day_separation, lower_CI,'LineWidth',2,'LineStyle','--','Color',[0 0 0]);
    hold off;
    maxr = .99; 
    tic_vect = [-maxr -.9 -.8 -.6 -.3 0 .3 .6 .8 .9 maxr];
    yticks(atanh(tic_vect));
    yticklabels(tic_vect);
    low_r = min(corr_mtx,[],'all');
    ylim(atanh([-0.3 maxr]));
    xlabel('Days Apart')
    ylabel('Rate Map Correlation (r)')
    title(title_string);
end

    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    

        

