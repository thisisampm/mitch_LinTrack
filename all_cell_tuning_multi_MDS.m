function [all_tcm, all_subj_idx, subj_cell] = all_cell_tuning_multi_MDS(subject_ids, folder_name, traces_type, reference_session, place_cell_only)
% plot the tuning curves of cells that are common across all session files
% sort rows by peak with reference to session number sort_by_sesh

% all_tcm = all_cell_tuning_multi([{'152-2'} {'159-2'}], 'control', reference_session)


% folder path
fp = fullfile('C:\Users\mitch\OneDrive - University of Toronto\MATLAB\mitch_LinTrack\', folder_name);

% details
num_spatial_bins = 40;
%traces_type = 2; % S, C, RAW

% Get the number of sessions from the number of files in the first folder
% besides the cellReg file

example_session_names = get_file_paths_all(fullfile(folder_name, subject_ids{1}));
cell_reg_index = contains(example_session_names, {'cell_regist','cellRegist'},'IgnoreCase',1);
sesh_nums = 1:numel(example_session_names(~cell_reg_index));

% iterate through subjects
%
subj_cell = cell(1,length(subject_ids));
subj_idx_cell = cell(1,length(subject_ids));
for isubj = 1:length(subject_ids)
    
    % get session file paths and cell reg matrix
    session_files = get_file_paths_targeted([fp '\' subject_ids{isubj}], {'.mat'});
    cell_regist_mtx_file = session_files(contains(session_files, {'cell_regist','cellRegist'}));
        variableInfo = who('-file',cell_regist_mtx_file{1});
        if ismember('cell_regist_mtx', variableInfo)
            load(cell_regist_mtx_file{1}, 'cell_regist_mtx')
        elseif ismember('cell_registered_struct', variableInfo)
            load(cell_regist_mtx_file{1}, 'cell_registered_struct')
            cell_regist_mtx = cell_registered_struct.cell_to_index_map;
        else
            error('bad cell_regist.mat file')
        end
    session_files = session_files(~contains(session_files, {'cell_regist_mtx','cellRegist'}));

    % compute tuning curves
    session_files(sesh_nums)
    [~, tuning_curve_matrix] = place_cell_tuning(session_files(sesh_nums), cell_regist_mtx(:,sesh_nums),traces_type, reference_session, place_cell_only);        
    
    % load
    subj_cell{isubj} = tuning_curve_matrix;
    subj_idx_cell{isubj} = repmat(isubj,size(tuning_curve_matrix,1),1);
    clearvars cell_regist_mtx tuning_curve_matrix;
end
%}
% save cell of all tuning curves
%save('subject_cell.mat', '-v7.3')
%load('C:\Users\ampm1\Documents\MATLAB\mitch_LinTrack\subject_cell_ctl2wk.mat')


% combine
all_tcm = [];
all_subj_idx = [];
for isesh = 1:length(subj_cell)
    all_tcm = [all_tcm; subj_cell{isesh}];
    all_subj_idx = [all_subj_idx; subj_idx_cell{isesh}];
    size(all_tcm);
end

% only include cells active during reference session
all_tcm = all_tcm(~isnan(all_tcm(:,1,reference_session)),:,:);

% sort by peak
[~,sort_idx] = sort_rows_by_peak(all_tcm(:,:,reference_session));
all_tcm = all_tcm(sort_idx,:,:);
all_subj_idx = all_subj_idx(sort_idx);

% plot
figure
for isesh = 1:length(sesh_nums)
   subplot(1,length(sesh_nums),isesh)
   imagesc(zscore_mtx(all_tcm(:,:,isesh)')')
   caxis([-1.5 1.5])
end
final_folder_name = regexp(folder_name,filesep,'split');
sgtitle(sprintf('%s All Cell Tuning Curve \n Created using CNMFE output: %g', final_folder_name{end},traces_type));


