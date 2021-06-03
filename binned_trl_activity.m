function [trl_act_mtx] = binned_trl_activity(behavior_mtx, trace_col_vect, num_bins, trials)
% comptues mean activity in each bin on each trial for a SINGLE neuron
%
% trace input is a column same size1 as behavior_mtx
% outputs mtx with trial rows and bin columns

% preallocate mtx
trl_act_mtx = nan(length(trials), num_bins);

%min max trial pos
min_max_xpos = [min(behavior_mtx(~isnan(behavior_mtx(:,4)),2))...
    max(behavior_mtx(~isnan(behavior_mtx(:,4)),2))];

% bins
bin_edges_prep = linspace(min_max_xpos(1), min_max_xpos(2), num_bins+1);
bin_edges = [bin_edges_prep(1:end-1)' bin_edges_prep(2:end)'];

% if vector of spikes, compute rate
if unique(trace_col_vect(trace_col_vect>0)) == 1 % spike vect test
    
    % iterate through trials
    for itrl = 1:length(trials)

        % trial index
        current_trl = trials(itrl);
        trl_idx = behavior_mtx(:,4)==current_trl;

        % iterate through bins
        for ibin = 1:num_bins
            bin_idx = behavior_mtx(:,2)>bin_edges(ibin,1) & behavior_mtx(:,2)<=bin_edges(ibin,2);
            trl_act_mtx(itrl,ibin) = sum(trace_col_vect(trl_idx & bin_idx))/(sum(trl_idx & bin_idx)/100); %hz
        end
    end

% else if vector of flor, compute mean
else
    % iterate through trials
    for itrl = 1:length(trials)

        % trial index
        current_trl = trials(itrl);
        trl_idx = behavior_mtx(:,4)==current_trl;

        % iterate through bins
        for ibin = 1:num_bins
            bin_idx = behavior_mtx(:,2)>bin_edges(ibin,1) & behavior_mtx(:,2)<=bin_edges(ibin,2);
            trl_act_mtx(itrl,ibin) = mean(trace_col_vect(trl_idx & bin_idx));
        end
    end
end

