%% Function that plots a transparent rectangle over the times when a trial is ongoing in an axis
% MDS 
% Inputs:
% 1. ax, the axis handle for the target
% 2. trial_markers, a unique trial x [start end] matrix


function [ax, lines] = plot_trial_patch(behavior_mtx, trace_mtx,varargin)

% Patch function takes two vectors (e.g. x and y) and plots a polygon
% between the vertices [x(n),y(n)] therefore need to create a matrix
% consisting of the paired vectors for each trial start and end paired to
% the min and max values of the target plot

[ax, lines] = plot_traces_offset(trace_mtx);
axis(ax,'padded');
[trials, trial_markers, dir_vect] = get_trials(behavior_mtx);

% preallocate
ntrial = length(trials);
x = nan(4, ntrial);
y = nan(4, ntrial);

% For each trial get the rectangle shape
for i = 1:ntrial
    x(:,i) = [trial_markers(i,1), trial_markers(i,1), trial_markers(i,2), trial_markers(i,2)]';
    y(:,i) = [ax.YLim(1) ax.YLim(2) ax.YLim(2) ax.YLim(1)]';
end

% Plot the polygons
rightward_trials = dir_vect == 0;
leftward_trials = dir_vect == 1;
patch(ax, x(:,rightward_trials), y(:,rightward_trials), 'r', 'FaceAlpha', 0.25,'LineStyle','none');
patch(ax, x(:,leftward_trials), y(:,leftward_trials), 'b', 'FaceAlpha', 0.25,'LineStyle','none');





