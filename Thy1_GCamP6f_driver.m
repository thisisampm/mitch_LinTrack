%%

%% Get files
%%
files_1881 = get_file_paths_targeted('exercise/188-1','2021-');
cellReg_1881 = get_file_paths_targeted('exercise/188-1','cellRegist');
load(cellReg_1881{1});
ind_1881 = cell_registered_struct.cell_to_index_map;
%%
files_200_4 = get_file_paths_targeted('exercise/200-4','2021-');
cellReg_200_4 = get_file_paths_targeted('exercise/200-4','cellRegist');
load(cellReg_200_4{1});
ind_200_4 = cell_registered_struct.cell_to_index_map;
%%
files_200_5 = get_file_paths_targeted('exercise/200-5','2021-');
cellReg_200_5 = get_file_paths_targeted('exercise/200-5','cellRegist');
load(cellReg_200_5{1});
ind_200_5 = cell_registered_struct.cell_to_index_map;
%%
files_196_3 = get_file_paths_targeted('control/196-3','2021-');
cellReg_196_3 = get_file_paths_targeted('control/196-3','cellRegist');
load(cellReg_196_3{1});
ind_196_3 = cell_registered_struct.cell_to_index_map;

%%
files_199_3 = get_file_paths_targeted('irradiation/199-3','2021-');
cellReg_199_3 = get_file_paths_targeted('irradiation/199-3','cellRegist');
load(cellReg_199_3{1});
ind_199_3 = cell_registered_struct.cell_to_index_map;

%% Run exercise animal 188-1
[tuning_mtx_ref_1881, tuning_mtx_all_1881] = all_cell_tuning(files_1881, ind_1881,1);
[cor_mtx_1881, ebp_reference_1881, ebp_adjacent_1881] = sesh_to_sesh_corrs(tuning_mtx_all_1881, 1);

%% Run control animal 196-3
[tuning_mtx_ref_196_3, tuning_mtx_all_196_3] = all_cell_tuning(files_196_3, ind_196_3,1);
[cor_mtx_196_3, ebp_reference_196_3, ebp_adjacent_196_3] = sesh_to_sesh_corrs(tuning_mtx_all_196_3, 1);

%% Run exercise animal 200-4
[tuning_mtx_ref_200_4, tuning_mtx_all_200_4] = all_cell_tuning(files_200_4, ind_200_4,1);
[cor_mtx_200_4, ebp_reference_200_4, ebp_adjacent_200_4] = sesh_to_sesh_corrs(tuning_mtx_all_200_4, 1);

%% Run exercise animal 200-5
[tuning_mtx_ref_200_5, tuning_mtx_all_200_5] = all_cell_tuning(files_200_5, ind_200_5,1);
[cor_mtx_200_5, ebp_reference_200_5, ebp_adjacent_200_5] = sesh_to_sesh_corrs(tuning_mtx_all_200_5, 1);

%% Run irradiation animal 199-3
[tuning_mtx_ref_199_3, tuning_mtx_all_199_3] = all_cell_tuning(files_199_3, ind_199_3,1);
[cor_mtx_199_3, ebp_reference_199_3, ebp_adjacent_199_3] = sesh_to_sesh_corrs(tuning_mtx_all_199_3, 1);

%% Run Irradiation Multi Animal plots
ctr_mice = [{'152-2'},{'159-2'},{'167-4'},{'196-3'},{'210-3'}];
irr_mice = [{'160-2'},{'167-2'},{'167-3'},{'190-3'},{'199-3'},{'209-2'},{'215-2'},{'215-3'}];
run_mice = [{'188-1'},{'200-2'},{'200-4'},{'200-5'}];

%% All cells spatial tuning
traces_type = 3;
reference_session = 1;
place_cell_only = 0;
[all_tcm_CTR, all_subj_idx_CTR, subj_cell_CTR] = all_cell_tuning_multi_MDS(ctr_mice, 'reload\control', traces_type, reference_session, place_cell_only);
[all_tcm_IRR, all_subj_idx_IRR, subj_cell_IRR] = all_cell_tuning_multi_MDS(irr_mice,'reload\irradiation',traces_type, reference_session, place_cell_only);
[all_tcm_RUN, all_subj_idx_RUN, subj_cell_RUN] = all_cell_tuning_multi_MDS(run_mice, 'reload\exercise', traces_type, reference_session, place_cell_only);

%%
reference_session = 1;
[cor_mtx_all_CTR, ebp_reference_all_CTR, ebp_adjacent_all_CTR] = sesh_to_sesh_corrs_multi_MDS(subj_cell_CTR, reference_session);
[cor_mtx_all_IRR, ebp_reference_all_IRR, ebp_adjacent_all_IRR] = sesh_to_sesh_corrs_multi_MDS(subj_cell_IRR, reference_session);
[cor_mtx_all_RUN, ebp_reference_all_RUN, ebp_adjacent_all_RUN] = sesh_to_sesh_corrs_multi_MDS(subj_cell_RUN, reference_session);

%% Place cell only spatial tuning
traces_type = 3;
reference_session = 1;
place_cell_only = 1;
[pc_all_tcm_CTR, pc_all_subj_idx_CTR, pc_subj_cell_CTR] = all_cell_tuning_multi_MDS([{'152-2'} {'159-2'} {'167-4'} {'196-3'} {'210-3'}], 'reload\control', traces_type, reference_session, place_cell_only);
[pc_all_tcm_IRR, pc_all_subj_idx_IRR, pc_subj_cell_IRR] = all_cell_tuning_multi_MDS([{'160-2'} {'167-2'} {'167-3'} {'190-3'} {'199-3'} {'215-2'} {'215-3'}],'reload\irradiation',traces_type, reference_session, place_cell_only);
[pc_all_tcm_RUN, pc_all_subj_idx_RUN, pc_subj_cell_RUN] = all_cell_tuning_multi_MDS([{'188-1'} {'200-2'} {'200-4'} {'200-5'} {'214-5'}], 'reload\exercise', traces_type, reference_session, place_cell_only);

%%
[cor_mtx_pc_CTR, ebp_reference_pc_CTR, ebp_adjacent_pc_CTR] = sesh_to_sesh_corrs_multi_MDS(pc_subj_cell_CTR, reference_session);
[cor_mtx_pc_IRR, ebp_reference_pc_IRR, ebp_adjacent_pc_IRR] = sesh_to_sesh_corrs_multi_MDS(pc_subj_cell_IRR, reference_session);
[cor_mtx_pc_RUN, ebp_reference_pc_RUN, ebp_adjacent_pc_RUN] = sesh_to_sesh_corrs_multi_MDS(pc_subj_cell_RUN, reference_session);

%% Population vector analysis 

groups = {'Control','Irradiation','Exercise'};

pv_ctr = pv_group_corr('reload\control');
pv_irr = pv_group_corr('reload\irradiation');
pv_ex = pv_group_corr('reload\exercise');

%% 

ctr1963 = get_file_paths_all('reload\control\196-3');
sc = [];
si = [];
for i = 1:numel(ctr1963(1:end-1))
    load(ctr1963{i});
    sc = [sc; stability_score(behavior_mtx,traces)];
    si = [si; information_score(behavior_mtx,traces,3)];
end

%% Days Apart Correlation Analysis
subj_corr_mtx_ctr = cell(size(ctr_mice));
for i = 1:numel(ctr_mice)
    subj_corr_mtx_ctr{i} = days_apart_corr(subj_cell_CTR{i});
end
subj_corr_mtx_irr = cell(size(irr_mice));
for i = 1:numel(irr_mice)
    subj_corr_mtx_irr{i} = days_apart_corr(subj_cell_IRR{i});
end
subj_corr_mtx_run = cell(size(run_mice));
for i = 1:numel(run_mice)
    subj_corr_mtx_run{i} = days_apart_corr(subj_cell_CTR{i});
end

all_corr_mtx_ctr = cat(1,subj_corr_mtx_ctr{:});
all_corr_mtx_irr = cat(1,subj_corr_mtx_irr{:});
all_corr_mtx_run = cat(1,subj_corr_mtx_run{:});
    
%% Temp

figure;
t = tiledlayout('flow');
for i = 1:length(weird_cells)
    nexttile(t)
    plot(all_tcm(weird_corr(i),:,1)); hold on; plot(all_tcm(weird_corr(i),:,2));hold off;
    title(string(weird_corr(i)));axis padded;
    plot_cell_activity_by_trial(rate_mtx_3d_day1,weird_cells(i,1));
    title(['Day 1 ' string(weird_corr(i))]);
    plot_cell_activity_by_trial(rate_mtx_3d_day2,weird_cells(i,2));
    title(['Day 2 ' string(weird_corr(i))]);
    plot_spike_positions(behavior_mtx_day1,traces_day1,weird_cells(i,1));
    title(['Day 1 ' string(weird_corr(i))]);
    plot_spike_positions(behavior_mtx_day2,traces_day2,weird_cells(i,2));
    title(['Day 2 ' string(weird_corr(i))]);
end




















































