%% Function to define the population vectors for each spatial bin during a session
function [pv_mtx, pv_mean] = pv_trial(behavior_mtx, traces, varargin)

if isempty(varargin)
    trace_type = 3;
    nbins = 40;
elseif length(varargin) == 1
    trace_type = varargin{1};
    nbins = varargin{1};
elseif length(varargin) == 2
    trace_type = varargin{1};
    nbins = varagin{2};
end

if trace_type == 1
    traces = traces(:,:,1);
    sampling_frequency = 20; % 20 Hz reload
elseif trace_type == 2
    traces = traces(:,:,2);
elseif trace_type == 3
    sampling_frequency = 20;
    z_threshold = 2;
    traces = extract_binary(traces(:,:,3),sampling_frequency,z_threshold);
end

% Cover continous position data to spatial bins along track

xbinned = discretize(behavior_mtx(:,2),40); % This iteration bins the 
% entire track, not just the portion of the track inlcuded in a trial,
% and is inconsistent with Adam's binned_trl_activty.m approach that bins
% the 80% of the track considered for a trial.  


% Get necessary indices
trials = unique(behavior_mtx(~isnan(behavior_mtx(:,4)),4)); % list of unique trials
bins = unique(xbinned(~isnan(xbinned))); % list of unique spatial bins

% Create a neuron x spatial bin x trial matrix to hold population vector
% data
n_neurons = size(traces,1);
ntrials = numel(trials);
pv_mtx = nan(n_neurons, nbins, ntrials);


%% For each trial, get activity of each neuron in each bin
for i = 1:numel(trials)
    trial_idx = behavior_mtx(:,4) == i; % Index for when mouse is running a given trial
    ca_trial = traces(:,trial_idx); % Ca flourescence data during trial
    pos_trial = xbinned(trial_idx); % Binned position data during trial
    for j = 1:numel(bins)
        in_bin = pos_trial == bins(j); % Index for time during the trial when in a given spatial bin
        ca_bin = ca_trial(:,in_bin); % Ca flourescence for each neuron in given spatial bin on that trial
        if trace_type == 1
            pv_mtx(:,j,i) = sum(ca_bin,2)./(sum(in_bin)/sampling_frequency);
        elseif trace_type == 2
            pv_mtx(:,j,i) = mean(ca_bin, 2,'omitnan'); % For each bin of each trial, create the PV as the vector of all neurons Ca flourescence at that bin
        % I just used flourescence, does that make sense?Or should be df/F?
        elseif trace_type == 3
            pv_mtx(:,j,i) = sum(ca_bin,2)/sum(in_bin);
        end
    end
end

%% Calculate the population vector that is the average of activity over the session
pv_mean = mean(pv_mtx, 3, 'omitnan');
    
    
