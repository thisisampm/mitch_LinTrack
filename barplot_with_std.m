%% Function that takes input data, plots the mean as bar and standard deviation as error.

function [bar_handle, err_handle] = barplot_with_std(data)

%% Calculate the mean and standard deviation of each group

if isnumeric(data) % If a matrix, consider columns to be groups
    group_means = mean(data,'omitnan'); % Calculate group means
    group_error = std(data, 'omitnan'); % Calculate group standard deviations
    ngroups = size(data,2); % Take the number of columns as the number of groups
elseif iscell(data)% If groups have a different number of values they can be inserted in a cell
    group_means = nan(1, numel(data)); % Preallocate a row vector to hold group means
    group_error = nan(1, numel(data)); % Preallocate row vector to hold standard deviation
    ngroups = numel(data); % The number of groups is a number of elements in the cell array
    for igroup = 1:numel(data)
        current_group_data = data{igroup}; % Pull the data out of the cell
        group_means(igroup) = mean(current_group_data, 'omitnan'); % Calculate mean
        group_error(igroup) = std(current_group_data, 'omitnan'); % Calculate standard deviation
    end
    
else
    error('First input is data in a subject x group matrix or a cell array of groups.')
end

%% Plot

xpos = 1:ngroups; % Create a position vector to align the groups along the x axis
figure; % Open a new figure
bar_handle = bar(xpos, group_means, 'FaceColor','flat'); % Create bar plot
hold on
err_handle = errorbar(xpos, group_means, group_error, group_error,'LineStyle','none'); % Create the errorbar plot
hold off

%% Aesthetics
group_colors = linspecer(ngroups); % Define distinct colours for each of the groups
bar_handle.CData = group_colors; % Set the bar graph colors

err_handle.CapSize = 10; % Increase the error bar horizontal width
err_handle.LineWidth = 2; % Make the errorbar bolder
err_handle.Color = [0 0 0]; % Force error bars to be black