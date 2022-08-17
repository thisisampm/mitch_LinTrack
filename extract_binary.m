function [binarized_trace, filtered_trace,norm_trace,d1_trace] = extract_binary(calcium_trace,sampling_frequency, z_threshold)
%extract_binary Converts raw calcium traces into binary traces
%   This function converts raw calcium traces inputs into a binarized output vector

% Copyright (C) 2017-2019 by Guillaume Etter
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or any
% later version.
% Contact: etterguillaume@gmail.com

%% Edited for use with a neuron x timestep matrix of neurons
% Mitch de Snoo - 12 May 2022

% calcium_trace: double/float matrix representing activity of a neurons.
% Each row represents 1 neuron and the columns are observations of Ca
% flourescence at that timestep
% 
% sampling frequency: frequency at which calcium imaging has been done (in Hz or fps)
%
% z_threshold: standard deviation threshold above which calcium activity is
% considered higher than noise

%% Parameters
[bFilt,aFilt] = butter(2,  2/(sampling_frequency/2), 'low');

filtered_trace = nan(size(calcium_trace)); % Preallocate 
ca_nan_pos = isnan(calcium_trace); % Identify NaNs created when first Ca timestamp is not necessarily equivalent to first behaviour timestamp

ncells = size(calcium_trace,1); % Number of neurons with distinct traces
traces_to_filter = reshape(calcium_trace(~ca_nan_pos),[ncells length(calcium_trace(~ca_nan_pos))/ncells]); % Reshape non-NaN traces into traces x timestamp 

filtered_trace(~ca_nan_pos) = (filtfilt(bFilt,aFilt,traces_to_filter'))'; % filtfilt operates along first dimension greater than 1, so need to transpose
norm_trace = filtered_trace./std(filtered_trace,[],2,'omitnan');
d1_trace = diff(filtered_trace,1,2); % Operate 1st derivative along each neuron
d1_trace(:,end+1) = 0;

binarized_trace = calcium_trace*0;
binarized_trace(norm_trace>z_threshold & d1_trace>0) = 1;
end

