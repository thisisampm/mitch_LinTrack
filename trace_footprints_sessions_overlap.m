function trace_footprints_sessions_overlap(footprint_cell, registration_mtx)
% plots outline of all footprints
% one color for first session, one for second, and one for overlapping
% cells (spatial footprints averaged)
%
% input is a cell matrix
% each cell contains 3d matrices of footprints from a single session
% neuron, ~length, ~width
%
% and the cell registration matrix output by cellreg
% cell_registered_struct.cell_to_index_map
% including only the columns corresponding the the cells in footprint_cell


%% place cells into new matrices
new_fp_cell = cell(1, length(footprint_cell)+1);
    
% not overall common
for icell = 1:length(footprint_cell)
    index_rows = registration_mtx(registration_mtx(:,icell)>0 & sum(registration_mtx==0,2)>0, icell);
    uncommon_footprint_cell{icell} = footprint_cell{icell}(index_rows,:,:);
end

% overall common
for icell = 1:length(footprint_cell)
    index_rows = registration_mtx(registration_mtx(:,icell)>0 & sum(registration_mtx==0,2)==0, icell);
    common_footprint_cell{icell} = footprint_cell{icell}(index_rows,:,:);
end



%% footprint figure
hold on

% plot colors
plot_colors = distinguishable_colors(length(footprint_cell));
plot_colors = [plot_colors; plot_colors./1.75];

% legend trick
legend_trick(plot_colors, 'o')

% plot
trace_footprints_sessions(uncommon_footprint_cell, plot_colors(1:length(uncommon_footprint_cell), :))
trace_footprints_sessions(common_footprint_cell, plot_colors(length(common_footprint_cell)+1 : end, :))

% legend
legend_string = [];
for iinput = 1:length(uncommon_footprint_cell)
    legend_string = [legend_string; ['Unique Session 0' num2str(iinput)]];
end
for iinput = 1:length(common_footprint_cell)
    legend_string = [legend_string; ['Common Session 0' num2str(iinput)]];
end
legend(legend_string)


%% Venn diagram
if ismember(length(footprint_cell), [2 3])

% session cell counts
session_cell_counts = nan(size(footprint_cell));
for icell = 1:length(footprint_cell)    
    session_cell_counts(icell) = sum(registration_mtx(:,icell)>0);
end

session_cell_counts

% overlapping cell count
common_cell_count = sum(sum(registration_mtx==0,2)==0)

% total cells
total_cells = sum(session_cell_counts)-common_cell_count

%figure
figure
[H,S] = venn(session_cell_counts, common_cell_count);
for ilabel = 1:size(S.ZoneCentroid,1)-1
    text(S.ZoneCentroid(ilabel,1), S.ZoneCentroid(ilabel,2), ...
        [legend_string(ilabel,:) ' cells (' num2str(session_cell_counts(ilabel)-common_cell_count) ';' num2str(round(((session_cell_counts(ilabel)-common_cell_count)/total_cells)*1000)/10) '%)'])
end
text(S.ZoneCentroid(end,1), S.ZoneCentroid(end,2), ['Common cells (' num2str(common_cell_count) ';' num2str(round((common_cell_count/total_cells)*1000)/10) '%)'])
axis equal off

end


