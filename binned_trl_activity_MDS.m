function [trl_act_mtx] = binned_trl_activity_MDS(behavior_mtx, trace_mtx, trace_type, nbins, trials)
% comptues mean activity in each bin on each trial for a neuron x time
% matrix
%
% Inputs: 
% 1. behavior_mtx
% 2. traces - formated as neurons x timestamps
% 3. Number of space bins
% 4. Index of which trials to include 
% Outputs:
% 1.mtx with trial rows and bin columns

sampling_frequency = 20; % Set to 20 Hz reload or 100 Hz for original AMPM loading

if trace_type == 1
    trace_mtx = trace_mtx(:,:,trace_type);
elseif trace_type == 2
    trace_mtx = trace_mtx(:,:,trace_type);
elseif trace_type == 3
    z_threshold = 2;
    trace_mtx = extract_binary(trace_mtx(:,:,trace_type),sampling_frequency,z_threshold);
end

%%
% preallocate mtx
ncells = size(trace_mtx,1);
ntrials = length(trials);
trl_act_mtx = nan(ntrials, nbins, ncells);

%min max trial pos
min_max_xpos = [min(behavior_mtx(~isnan(behavior_mtx(:,4)),2))...
    max(behavior_mtx(~isnan(behavior_mtx(:,4)),2))];

% bins
bin_edges_prep = linspace(min_max_xpos(1), min_max_xpos(2), nbins+1);
bin_edges = [bin_edges_prep(1:end-1)' bin_edges_prep(2:end)'];

%xpos_binned = discretize(behavior_mtx(:,4),nbins);

% if vector of spikes, compute rate
%if unique(trace_mtx(trace_mtx>0)) == 1 % spike vect test


for itrl = 1:ntrials
    % Identify which timesteps correspond to current trial
    current_trl = trials(itrl);
    trl_idx = behavior_mtx(:,4) == current_trl;
    % Iterate through bins calculating calcium activity within each
    for ibin = 1:nbins
        bin_idx = behavior_mtx(:,2)>bin_edges(ibin,1) & behavior_mtx(:,2)<=bin_edges(ibin,2);
        if trace_type == 1 % If spikes (s), compute rate
            trl_act_mtx(itrl,ibin,:) = sum(trace_mtx(:,trl_idx & bin_idx),2)./(sum(trl_idx & bin_idx)/sampling_frequency); %hz
        elseif trace_type == 2 % If C compute mean flouresence
            trl_act_mtx(itrl,ibin,:) = mean(trace_mtx(:,trl_idx & bin_idx),2,'omitnan');
        elseif trace_type == 3 % if C_raw, expect it to be binarized and compute fraction of active frames
            trl_act_mtx(itrl,ibin,:) = sum(trace_mtx(:,trl_idx & bin_idx),2)./(sum(trl_idx & bin_idx));            
        end
    end
end

% Indicate what trace type was used
%{
if trace_type == 1
    fprintf('Trial activity calculated using S and sampling rate of %g Hz. \n',sampling_frequency);
elseif trace_type == 2
    fprintf('Trial activity calculated using mean C flouresence. \n');
elseif trace_type == 3
    fprintf('Trial activity calculate using binarized C_raw. \n');
end
 %}           
%% Old AMPM code, can find in original binned trial activity 
%{
if trace_type == 1   
    % iterate through trials
    for itrl = 1:length(trials)

        % trial index
        current_trl = trials(itrl);
        trl_idx = behavior_mtx(:,4)==current_trl;

        % iterate through bins
        for ibin = 1:nbins
            bin_idx = behavior_mtx(:,2)>bin_edges(ibin,1) & behavior_mtx(:,2)<=bin_edges(ibin,2);
            trl_act_mtx(itrl,ibin,:) = sum(trace_mtx(:,trl_idx & bin_idx),2)./(sum(trl_idx & bin_idx)/100); %hz
        end
    end

% else if vector of flor, compute mean
elseif trace_type == 2
    % iterate through trials
    for itrl = 1:length(trials)

        % trial index
        current_trl = trials(itrl);
        trl_idx = behavior_mtx(:,4)==current_trl;
        
        % Count bins vectorized
        x_binned = discretize(behavior_mtx(trl_idx,2),40); % Bin the position data during trial
        %x_binned = 
        %trl_act_mtx(itrl,:) = 
        

        % iterate through bins
        for ibin = 1:nbins
            bin_idx = behavior_mtx(:,2)>bin_edges(ibin,1) & behavior_mtx(:,2)<=bin_edges(ibin,2);
            trl_act_mtx(itrl,ibin,:) = mean(trace_mtx(:,trl_idx & bin_idx),2,'omitnan');
        end        
    end        
end
%}
