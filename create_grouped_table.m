function [grouped_data_table] = create_grouped_table(data, varargin)

% Function that takes cell tuning correlation matrices and reshapes them 

%%%%%%%% Inputs %%%%%%%%%%%
% Data: a cell array with the corr_mtx for each group output from
% sesh_to_sesh_corrs_multi
% OPTIONAL:
% group_names: a 1D cell array of names corresponding to the names of
% the groups in data
% mouse_id: group x a cell array

ngroups = length(data);
% If not given group names, categorize as numbers.
if ~isempty(varargin)
    if length(varargin) == 1
        group_names = categorical(varargin{1});
    elseif length(varargin) == 2
        group_names = categorical(varargin{1});
        mouse_id_cell = varargin{2};
    else
        error('Too many inputs')
    end
else
    group_names = categorical(1:ngroups);
end
    

%% 
for igroup = 1:ngroups
    corr_mtx = data{igroup}; % Extract day x day x animal matrix of average correlations
    corr_with_day1 = squeeze(corr_mtx(1,:,:)); % Restrict to top row for each animal, which represents correlations with Day 1
    nanimals = size(corr_with_day1,2); % Each column represents 1 animal
    nsessions = size(corr_with_day1,1); % Each row represents average spatial correlation with day 1.
    corr_vector = corr_with_day1(:); % Reshape corr_mtx into a column vector
    group_vector = repelem(group_names(igroup),numel(corr_vector),1); % Categorical column vector indicating which group each corr sample belongs to
    session_vector = (repmat((1:2:nsessions*2-1)',nanimals,1)); % Column vector indicating which day n is correlated with day 1
    % If not given mouse ids, categorize them as mouse 1 - mouse nmice    
    if exist('mouse_id_cell','var')
        mouse_ids = mouse_id_cell{igroup};
    else
        mouse_ids = categorical(1:nanimals);
    end
    animal_vector = repmat(mouse_ids, nsessions, 1); % Categorical vector indicating which mouse each corr data point belongs to
    animal_vector = animal_vector(:);
    igroup_table = table(corr_vector,group_vector,session_vector,animal_vector,...
        'VariableNames',{'Corr','Group','Day','Mouse'}); % Create table for current group
    if igroup == 1
        grouped_data_table = igroup_table;
    else
        grouped_data_table = [grouped_data_table; igroup_table];% Append the current group to the overall table
    end
end

    
    
    
    
    
    
    
