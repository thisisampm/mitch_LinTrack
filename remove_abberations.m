function traces = remove_abberations(traces, max_change)
% removes isolated values from C_raw that are more different than their
% surrounding pixels by at least max_change. Then interps a new value.

%figure; hold on; plot(traces)
traces_z = ( traces- mean(traces) ) ./std(traces);

% absolute changes
chng = abs(traces_z(2:end)-traces_z(1:end-1));

% find pairs of values greater than max_change
gtm = chng>max_change;
gtm = (gtm(1:end-1) + gtm(2:end))==2;

% bad pixels
bp = find(gtm==1)+1;

% interpolate
lt = 1:length(traces);
traces(bp) = nan; % delete from raw
nnan_idx = ~isnan(traces);
traces = interp1(lt(nnan_idx),traces(nnan_idx),lt);
