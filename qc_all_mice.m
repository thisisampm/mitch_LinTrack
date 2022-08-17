% Script to generate a quality control table

groups = {'control','irradiation','exercise'};

for igroup = 1:numel(groups)
    [mouse_folders, mice] = get_folder_paths_all(groups{igroup});
    nmice = numel(mouse_folders);
    for imouse = 1:nmice
        disp(['Computing quality control metrics for ', groups{igroup},' ', mice{imouse}])
        mouse_qc_table = quality_control(mouse_folders{imouse});
        if igroup == 1 && imouse == 1
            qc_table = mouse_qc_table;
        else
            qc_table = [qc_table; mouse_qc_table];
        end
    end
end
            
        