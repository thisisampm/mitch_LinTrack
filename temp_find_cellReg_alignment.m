%% Check that the away that cell reg is used in all_cell_tuning.m doesn't result in mismatched cells.

row424 = cell_registered_struct.cell_to_index_map(424,:); % row that represents session1, cell172, and it's matched cell in later sessions
row172 = cell_registered_struct.cell_to_index_map(172,:); % 