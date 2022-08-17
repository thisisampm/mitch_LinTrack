
reference_session = 1;

%% Get all relevant mice for each experimental group.

group_names = [{'control'},{'irradiation'},{'exercise'}];

to_exclude = {'G7','old','Old'};
%{
[ctr_folder_paths, ctr_folder_names] = get_folder_paths_all('control',1);
ctr_folder_paths = ctr_folder_paths(~contains(ctr_folder_paths,to_exclude));
ctr_folder_names = ctr_folder_names(~contains(ctr_folder_names,to_exclude));

[irr_folder_paths, irr_folder_names] = get_folder_paths_all('irradiation',1);
irr_folder_paths = irr_folder_paths(~contains(irr_folder_paths,to_exclude));
irr_folder_names = irr_folder_names(~contains(irr_folder_names,to_exclude));

[run_folder_paths, run_folder_names] = get_folder_paths_all('exercise',1);
run_folder_paths = run_folder_paths(~contains(run_folder_paths,to_exclude));
run_folder_names = run_folder_names(~contains(run_folder_names,to_exclude));
%}
%% Calculate information content, place cell, and place field peak
% Col 1 Group index (control = 1, irradiation = 2, exercise = 3)
% Col 2 Mouse index (number within group i.e. 1 for control 152-2 and
% irradiation 160-2
% Col 3 Session index
% Col 4 Spatial information content
% Col 5 Place cell criteria met, logical index
% Col 6 Centroid of place field peak

if ~exist('place_cell_info','var')

    place_cell_info = [];


    for igroup = 1:numel(group_names)
        [folder_paths, folder_names] = get_folder_paths_all(group_names{igroup},1);
        folder_paths = folder_paths(~contains(folder_paths,to_exclude));
        folder_names = folder_names(~contains(folder_names,to_exclude));
        for imouse = 1:numel(folder_names)
            all_files = get_file_paths_all(folder_paths{imouse});
            session_files = all_files(~contains(all_files,{'cell_regist','cellRegis'}));
            for isession = 1:numel(session_files)
                tic
                load(session_files{isession});
                if ~exist('place_cell_mtx','var')
                    [IC, pc_idx, place_field_peak] = information_score(behavior_mtx,traces);
                    disp(['Finished calculating ', session_files{isession}])
                elseif exist('place_cell_mtx','var')
                    IC = place_cell_mtx(:,1);
                    pc_idx = place_cell_mtx(:,2);
                    place_field_peak = place_cell_mtx(:,3);
                    disp(['Loading place cell data from ' session_files{isession}])
                end
                group_index = repmat(igroup,length(IC),1); % Column vector indicating group number
                mouse_index = repmat(imouse,length(IC),1); % Column vector indicating mouse number
                session_index = repmat(isession,length(IC),1); % Column vector indicating session number
                session_output = [group_index, mouse_index, session_index, IC, pc_idx, place_field_peak];
                place_cell_info = [place_cell_info; session_output];
                
                toc
            end
        end
    end
    
end

%% Calculate the fraction of place cells for each group

ctr_idx = place_cell_info(:,1) == 1;
irr_idx = place_cell_info(:,1) == 2;
run_idx = place_cell_info(:,1) == 3;

n_all_cells = nan(1,numel(group_names)); % 
n_pc = nan(1,numel(group_names)); % 

bins = unique(place_cell_info(:,6)); % Get the spatial bins while running trials
all_cell_peak_distribution = nan(numel(group_names),length(bins)); % Preallocate a group x bin number for all cell peak activity
pc_peak_distribution = nan(numel(group_names),length(bins)); % Preallocate a group x bin matrix for place cell peak activity

 
for i = 1:numel(group_names)
    group_idx = place_cell_info(:,1) == i; % Index for cells belonging to each of the three groups
    group_place_cell_info = place_cell_info(group_idx,:); % Select the place cell info for that group
    n_all_cells(i) = size(group_place_cell_info,1); % Number of columns represents every cell
    n_pc(i) = sum(group_place_cell_info(:,5)); % Sum the logical index indicating place cell to get total number
    
    all_cell_peak_distribution(i,:) = sum(group_place_cell_info(:,6) == bins'); % Find the number of peak fields for each bin for all cells
    group_pc_idx = logical(group_place_cell_info(:,5));
    pc_peak_distribution(i,:) = sum(group_place_cell_info(group_pc_idx,6) == bins'); % Find the number of peak fields for each bin for all cells
    
end

fraction_place_cells = n_pc./n_all_cells

%% Plot the distribution of all cell peak spatial firing
hex_colors = [{'#010101','#ED1C24','#3A54A5'}];
rgb_colors = hex2rgb(hex_colors);
close all;
figure;
subplot(2,1,1);
n_cells = sum(all_cell_peak_distribution,2);
all_cells_peak_fraction = 100*all_cell_peak_distribution./n_cells;
%{
bar(all_cells_peak_fraction(1,:));
hold on
bar(all_cells_peak_fraction(2,:));
bar(all_cells_peak_fraction(3,:));

%}0
bar(all_cells_peak_fraction','BarWidth',1.5)
title('Distribution of Peak Cell Firing for All Cells')
ylabel('% of Spatial Peaks')
xlabel('Spatial Bin')
legend(group_names);


%% Plot the distribution of all place cell peak spatial firing

n_pc = sum(pc_peak_distribution,2);
pc_peak_fraction = 100*pc_peak_distribution./n_pc;
%{
bar(pc_peak_fraction(1,:));
hold on
bar(pc_peak_fraction(2,:));
bar(pc_peak_fraction(3,:));
legend(group_names);
%}
subplot(2,1,2)
bar(pc_peak_fraction','BarWidth',1.5)
title('Distribution of Peak Cell Firing for Place Cells')
ylabel('% of Spatial Peaks')
xlabel('Spatial Bin')
legend(group_names);



%% Calculate the information content of all cells and all place cells for each group

bar_colors = flip([ 0.9290 0.6940    0.1250; 0.8500    0.3250    0.0980; 0    0.4470    0.7410]);


all_cell_ic = cell(1,numel(group_names));
all_cell_ic_mean = nan(1,numel(group_names));
all_cell_ic_std = nan(1,numel(group_names));

pc_ic = cell(1,numel(group_names));
pc_ic_mean = nan(1,numel(group_names));
pc_ic_std = nan(1,numel(group_names));



for igroup = 1:numel(group_names)
    group_index = place_cell_info(:,1) == igroup;
    all_cell_ic{igroup} = place_cell_info(group_index,4);
    all_cell_ic_mean(igroup) = mean(all_cell_ic{igroup},'omitnan');
    all_cell_ic_std(igroup) = std(all_cell_ic{igroup},'omitnan');
    
    group_pc_idx = logical(place_cell_info(group_index,5) == 1);
    pc_ic{igroup} = place_cell_info(group_pc_idx,4);
    pc_ic_mean(igroup) = mean(pc_ic{igroup},'omitnan');
    pc_ic_std(igroup) = std(pc_ic{igroup},'omitnan');
end

figure;
subplot(2,1,1)
b1 = bar(100*fraction_place_cells,'FaceColor','flat');
b1.CData = bar_colors;
xticklabels(group_names)
ylabel('% Place Cells')
axis([0 4 0 70])

subplot(2,1,2)
b2 = bar(pc_ic_mean,'FaceColor','flat');
b2.CData = bar_colors;
xticklabels(group_names);
ylabel(['Place Cell', newline,'Spatial Information (bits)'])
hold on
err = errorbar(1:3,pc_ic_mean, pc_ic_std,pc_ic_std);
err.Color = [0 0 0];
err.LineStyle = 'none';
hold off
}%