%% Trial number

% Load trial number data
nTrials_Control = [get_nTrials('control\152-2'),...
    get_nTrials('control\159-2'),...
    get_nTrials('control\167-4')];

nTrials_irradation = [get_nTrials('irradiation\160-2'),...
    get_nTrials('irradiation\167-2'),...
    get_nTrials('irradiation\167-3')];

% Plot as average trial number completed per session
bar_ax1 = make_bar_plot(nTrials_Control, nTrials_irradation);
bar_ax1.YLabel.String = 'Average Number of Trials Completed';

% Create line plot of trials each session
line_ax1 = make_line_plot(nTrials_Control, nTrials_irradation);
line_ax1.XTickLabel = cellstr(string(1:2:15));
line_ax1.YLabel.String = 'Number of Trials Completed';
line_ax1.XLabel.String = 'Day';
%% Trial duration

% Load trial duration data
trial_durations_control = [get_trial_duration('control\152-2'),...
    get_trial_duration('control\159-2'),...
    get_trial_duration('control\167-4')];

trial_durations_irradation = [get_trial_duration('irradiation\160-2'),...
    get_trial_duration('irradiation\167-2'),...
    get_trial_duration('irradiation\167-3')];

% Plot as average trial number completed per session
bar_ax2 = make_bar_plot(trial_durations_control, trial_durations_irradation);
bar_ax2.YLabel.String = 'Average Trial Duration (s)';

% Create line plot of trials each session
line_ax2 = make_line_plot(trial_durations_control, trial_durations_irradation);
line_ax2.XTickLabel = cellstr(string(1:2:15));
line_ax2.YLabel.String = 'Trial Duration (s)';
line_ax2.XLabel.String = 'Day';

%% Functions

function [nTrials] = get_nTrials(mouse_dir)
    sesh_info = dir(mouse_dir); % get info for each session 
    mask = ismember({sesh_info.name},{'.','..','cell_regist.mat'}); % exclude . and .. directories and cellReg file
    sesh_info(mask) = [];
    nTrials = nan(size(sesh_info)); % Preallocate matrix to hold number of trials for each session
    for i = 1:size(sesh_info,1) % For each session load behavior matrix and get the highest trial number (column 4)
        load(fullfile(sesh_info(i).folder,sesh_info(i).name), 'behavior_mtx')
        nTrials(i) = max(behavior_mtx(:,4));
    end
end

function [vel_by_trial] = get_velocity(mouse_dir)
    sesh_info = dir(mouse_dir); % get info for each session 
    mask = ismember({sesh_info.name},{'.','..','cell_regist.mat'}); % exclude . and .. directories and cellReg file
    sesh_info(mask) = [];
    vel_by_trial = nan(size(sesh_info)); % Preallocate matrix to hold average velocity during trials for each session
    for i = 1:size(sesh_info,1)
        load(fullfile(sesh_info(i).folder,sesh_info(i).name), 'behavior_mtx')
        intrial = ~isnan(behavior_mtx(:,4)); % determine when mouse is running a trial
        vel_by_trial(i) = mean(behavior_mtx(intrial,5),'omitnan');
    end
end

function [avg_trial_dur] = get_trial_duration(mouse_dir)
    sesh_info = dir(mouse_dir); % get info for each session 
    mask = ismember({sesh_info.name},{'.','..','cell_regist.mat'}); % exclude . and .. directories and cellReg file
    sesh_info(mask) = [];
    avg_trial_dur = nan(size(sesh_info)); % Preallocate matrix to hold average velocity during trials for each session
    for i = 1:size(sesh_info,1)
        load(fullfile(sesh_info(i).folder,sesh_info(i).name), 'behavior_mtx')
        trials = unique(behavior_mtx(~isnan(behavior_mtx(:,4)),4));
        trial_dur = nan(size(trials));        
        for j = 1:length(trials)
            %if ismember(j,unique(behavior_mtx(~isnan(behavior_mtx(:,4)),4)))
                trial_start = behavior_mtx(find(behavior_mtx(:,4)==trials(j),1,'first'),1);
                trial_end = behavior_mtx(find(behavior_mtx(:,4)==trials(j),1,'last'),1);
                trial_dur(j) = trial_end - trial_start;
            %end
        end
        avg_trial_dur(i) = mean(trial_dur);
    end
end

%%
function [bar_ax] = make_bar_plot(ctr_data, irr_data)
    figure;
    bar([mean(mean(ctr_data)),mean(mean(irr_data))],'LineWidth',2,'FaceAlpha',0); % bar plot means
    bar_ax = gca;
    hold on
    scatter(1,mean(ctr_data),'LineWidth',1.5,'jitter','on','jitterAmount',0.2) % Add individual control animals
    scatter(2,mean(irr_data),'LineWidth',1.5,'jitter','on','jitterAmount',0.2) % Add individual irradiation animals
    err_ctr = std(mean(ctr_data));
    errbar_control = errorbar(1, mean(mean(ctr_data)),err_ctr, err_ctr); % Add control std as errorbar
    errbar_control.Color = [0 0 0];
    errbar_control.LineStyle = 'none';
    errbar_control.LineWidth = 1.5;
    err_irr = std(mean(irr_data));
    errbar_irr = errorbar(2, mean(mean(irr_data)),err_irr, err_irr); % Add irradiation std as errorbar
    errbar_irr.Color = [0 0 0];
    errbar_irr.LineStyle = 'none';
    errbar_irr.LineWidth = 1.5;
    hold off
    ax1 = gca;
    ax1.XTickLabel = {'Control', 'Irradation'};
end

function [ax_line] = make_line_plot(ctr_data,irr_data)
    figure;
    ctr_line = plot(mean(ctr_data,2),'k','LineWidth',2); % Plot mean trial number over sessions for control
    axis([-inf inf 0 inf])
    hold on
    err_ctr = std(ctr_data,0,2);
    errbar_control = errorbar(1:size(ctr_data,1), mean(ctr_data,2),err_ctr, err_ctr); % Add control std as errorbar
    errbar_control.Color = [0 0 0];
    errbar_control.LineStyle = 'none';
    errbar_control.LineWidth = 1;
    plot(ctr_data,'k','LineStyle','--')

    irr_line = plot(mean(irr_data,2),'r','LineWidth',2); % Plot mean trial number over sessions for control
    err_irr = std(irr_data,0,2);
    errbar_control = errorbar(1:size(irr_data,1), mean(irr_data,2),err_irr, err_irr); % Add control std as errorbar
    errbar_control.Color = 'r';
    errbar_control.LineStyle = 'none';
    errbar_control.LineWidth = 1;
    plot(irr_data,'r','LineStyle','--')
    axis padded
    ax_line = gca;
    legend([ctr_line, irr_line],{'Control','Irradiated'});
    ax_line.Legend.Position = [0.1732 0.7877 0.1821 0.0821];
    hold off
end








