function all_tcm = all_cell_tuning_multi(subject_ids, folder_name, reference_session)
% plot the tuning curves of cells that are common across all session files
% sort rows by peak with reference to session number sort_by_sesh

% all_tcm = all_cell_tuning_multi([{'152-2'} {'159-2'}], 'control', reference_session)


% folder path
fp = ['C:\Users\ampm1\Documents\MATLAB\mitch_LinTrack\data\' folder_name];

% details
num_spatial_bins = 40;
traces_type = 3; % S, C, RAW
sesh_nums = 1:8;


% iterate through subjects
%{
subj_cell = cell(1,length(subject_ids));
for isubj = 1:length(subject_ids)
    
    % get session file paths and cell reg matrix
    session_files = get_file_paths_targeted([fp '\' subject_ids{isubj}], {'.mat'}); 
    cell_regist_mtx_file = session_files(contains(session_files, 'cell_regist'));
        load(cell_regist_mtx_file{1}, 'cell_regist_mtx')
    session_files = session_files(~contains(session_files, 'cell_regist_mtx'));

    % compute tuning curves
    session_files(sesh_nums)
    [~, tuning_curve_matrix] = all_cell_tuning(session_files(sesh_nums), cell_regist_mtx(:,sesh_nums), 1);
    
    % load
    subj_cell{isubj} = tuning_curve_matrix;
    clearvars cell_regist_mtx;
end
%}
% save cell of all tuning curves
%save('subject_cell.mat', '-v7.3')
load('C:\Users\ampm1\Documents\MATLAB\mitch_LinTrack\subject_cell_ctl2wk.mat')


% combine
all_tcm = [];
for isesh = 1:length(subj_cell)
    all_tcm = [all_tcm; subj_cell{isesh}];
    size(all_tcm)
end

% only include cells active during reference session
all_tcm = all_tcm(~isnan(all_tcm(:,1,reference_session)),:,:);

% sort by peak
[~,sort_idx] = sort_rows_by_peak(all_tcm(:,:,reference_session));
all_tcm = all_tcm(sort_idx,:,:);

% plot
figure
for isesh = 1:length(sesh_nums)
   subplot(1,length(sesh_nums),isesh)
   imagesc(zscore_mtx(all_tcm(:,:,isesh)')')
   caxis([-1.5 1.5])
end


