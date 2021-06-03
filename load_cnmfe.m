function [traces, frame_times] = load_cnmfe(fp)
% load flor traces and timestamps from cnmfe output

% load traces
fp

cnmfe_file = get_file_paths_targeted(fp, {'cnmfe'})
warning('off'); load(cnmfe_file{1}, 'st'); warning('on'); % andrew did something annoying

% set 'spike' output to binary
st.S(st.S>0) = 1;

% assign to traces output
traces = nan(size(st.C_raw,1), size(st.C_raw,2), 3);
traces(:,:,1) = st.S;
C_hold = st.C; C_hold(C_hold<=0.001)=0; traces(:,:,2) = C_hold;
traces(:,:,3) = zscore_mtx((st.C_raw)')';

% correct curve of raw traces
for itrace = 1:size(traces(:,:,3), 1)
   traces(itrace,:,3) = traces_mean_correct(traces(itrace,:,3), 300);
end

% load frame times
pkl_file = get_file_paths_targeted(fp, 'Ca_', 'ime');
frame_times = readmatrix(pkl_file{1});










