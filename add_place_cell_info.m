%% Function to add a place cell info matrix variable to the .mat files for each linear track day

function [] = add_place_cell_info(mouse_folder,trace_type)

%% Get all of the files to analyze place cell info
file_names = get_file_paths_all(mouse_folder); % Get all file names
cell_reg_file_index = contains(file_names, {'cellReg','Cell_Reg'},'IgnoreCase',1); % Find the cell reg file
file_names = file_names(~cell_reg_file_index); % Exclude cell reg file

%% Calculate the place cell info and save it to the .mat file

for isesh = 1:numel(file_names)
    clearvars behavior_mtx traces place_cell_mtx % Clear variables from previous iteration
    warning('off'); load(file_names{isesh},'behavior_mtx','traces','place_cell_mtx'); warning('on');
    if exist('place_cell_mtx','var') % Skip to next file if already contains a place cell mtx output
        warning([file_names{isesh} ,' already contains a place cell info variable.'])
        continue
    else
        disp(['Calculating place cell statistics for ', file_names{isesh}])
        [IC, pc_idx, place_field_peak] = information_score(behavior_mtx,traces,trace_type); % Calculate information content, index indicating if place cell criteria met, and place field centroid for each cell
        place_cell_mtx = [IC pc_idx place_field_peak]; % Concatenate into a matrix for saving as the output 
        % Col 1: information content
        % Col 2: Place cell index, 0 = Place cell 1 = Place cell
        % Col 3: Place field centroid
        save(file_names{isesh},'behavior_mtx','traces','place_cell_mtx')
    end    
end
        
        

        
