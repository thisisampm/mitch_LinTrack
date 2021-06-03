function trace_footprints_sessions(footprint_cell, varargin)
% plots outline of all footprints
% footprints are a cell matrix
% each cell contains 3d matrices of footprints from a single session
% neuron, ~length, ~width

% figure
hold on

% input
if ~isempty(varargin)
    plot_colors = varargin{1};
else
    plot_colors = distinguishable_colors(length(footprint_cell));
end

% legend trick
legend_trick(plot_colors, 'o')

% iterate through each cell
for icell = 1:length(footprint_cell)
    trace_footprints(footprint_cell{icell}, plot_colors(icell,:));
end

% aesthetics
axis([200 600 75 400]);

% legend
legend_string = [];
for iinput = 1:length(footprint_cell)
    legend_string = [legend_string; ['Session 0' num2str(iinput)]];
end
legend(legend_string)



