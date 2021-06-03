function traces_c = traces_mean_correct(trace_uc, smooth_window)
% corrects slope of trace by flattening a grossly smoothed trace

% remove aberrations
for i = 1:size(trace_uc,1)
    trace_uc(i,:) = remove_abberations(trace_uc(i,:), 3);
end

% smoothed trace
st = smooth(trace_uc, smooth_window)';

% mean of smoothed trace
mt = nanmean(st);

% difference
dt = st-mt;

% correct
traces_c = trace_uc-dt;

% norm
traces_c = (traces_c-min(traces_c))./max(traces_c-min(traces_c));
