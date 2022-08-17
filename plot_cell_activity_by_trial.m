% Function that plots the trial by trial activity of each cell selected
% Input:
% rate_mtx_3d: a trial x bin x cell activity matrix

function plot_cell_activity_by_trial(rate_mtx_3d,varargin)

% Take an optional index of cells to plot
if ~isempty(varargin)
    cell_ids = varargin{1};
    ncells = length(cell_ids);
else
    ncells = size(rate_mtx_3d,3);
    cell_ids = 1:ncells;
end


f = figure;
t = tiledlayout(f,'flow');
t.TileSpacing = 'tight';
t.XLabel.String = 'Spatial Bins';
t.YLabel.String = 'Trial Number';

for i = 1:ncells
    nexttile
    imagesc(rate_mtx_3d(:,:,cell_ids(i)));
    title(sprintf('Cell %g',cell_ids(i)));
end
    