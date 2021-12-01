function [all_tcm, all_subj_idx, subj_cell] = all_cell_tuning_multi(subject_ids, folder_name, reference_session)
% plot the tuning curves of cells that are common across all session files
% sort rows by peak with reference to session number sort_by_sesh

% all_tcm = all_cell_tuning_multi([{'152-2'} {'159-2'}], 'control', reference_session)


% folder path
fp = ['C:\Users\ampm1\Documents\MATLAB\mitch_LinTrack\data\' folder_name];

% details
num_spatial_bins = 40;
traces_type = 3; % S, C, RAW

if contains(folder_name, '3wk')
    sesh_nums = 1:12;
elseif contains(folder_name, 'extraday')
    sesh_nums = 1:9;
else
    sesh_nums = 1:8;
end

% iterate through subjects
%
subj_cell = cell(1,length(subject_ids));
subj_idx_cell = cell(1,length(subject_ids));
for isubj = 1:length(subject_ids)
    
    % get session file paths and cell reg matrix
    session_files = get_file_paths_targeted([fp '\' subject_ids{isubj}], {'.mat'});
    cell_regist_mtx_file = session_files(contains(session_files, 'cell_regist'));
        variableInfo = who('-file',cell_regist_mtx_file{1});
        if ismember('cell_regist_mtx', variableInfo)
            load(cell_regist_mtx_file{1}, 'cell_regist_mtx')
        elseif ismember('cell_registered_struct', variableInfo)
            load(cell_regist_mtx_file{1}, 'cell_registered_struct')
            cell_regist_mtx = cell_registered_struct.cell_to_index_map;
        else
            error('bad cell_regist.mat file')
        end
    session_files = session_files(~contains(session_files, 'cell_regist_mtx'));

    % compute tuning curves
    session_files(sesh_nums)
    [~, tuning_curve_matrix] = all_cell_tuning(session_files(sesh_nums), cell_regist_mtx(:,sesh_nums), 1);
    
    % load
    subj_cell{isubj} = tuning_curve_matrix;
    subj_idx_cell{isubj} = repmat(isubj,size(tuning_curve_matrix,1),1);
    clearvars cell_regist_mtx tuning_curve_matrix;
    
    %subj_idx_cell{isubj}
    
end
%}
% save cell of all tuning curves
% save('subject_cell.mat', '-v7.3')
% load('C:\Users\ampm1\Documents\MATLAB\mitch_LinTrack\subject_cell_ctl2wk.mat')


% combine
all_tcm = [];
all_subj_idx = [];
for isubj = 1:length(subject_ids)
    all_tcm = [all_tcm; subj_cell{isubj}];
    all_subj_idx = [all_subj_idx; subj_idx_cell{isubj}];
end



% only include cells active during reference session
all_tcm = all_tcm(~isnan(all_tcm(:,1,reference_session)),:,:);
all_subj_idx = all_subj_idx(~isnan(all_tcm(:,1,reference_session)));

% sort by peak
[~,sort_idx] = sort_rows_by_peak(all_tcm(:,:,reference_session));
all_tcm = all_tcm(sort_idx,:,:);

size(all_subj_idx)
all_subj_idx = all_subj_idx(sort_idx);
size(all_subj_idx)

% plot
figure
for isesh = 1:length(sesh_nums)
   subplot(1,length(sesh_nums),isesh)
   imagesc(zscore_mtx(all_tcm(:,:,isesh)')')
   caxis([-1.5 1.5])
end


