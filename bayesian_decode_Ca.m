function [class, posterior, tuning_curves] = bayesian_decode_Ca(sample, training, group, varargin)
% Bayesian decoder for spatial position based on Etter, Manseau &
% Williams (2020) "A Probabilistic Framework for Decoding Behavior From in vivo 
% Calcium Imaging Data"

%                                              Mitch de Snoo, 14 Apr 2022

% [CLASS, POSTERIOR] = BAYESIAN_DECODE(SAMPLE,TRAINING,GROUP, BINS) 
% classifies each row of the data in SAMPLE into one of the groups in TRAINING.  
% SAMPLE and TRAINING must be matrices with the same number of columns.  
% GROUP is a grouping variable for TRAINING.  Its unique values define 
% groups, and each element defines which group the corresponding row of 
% TRAINING belongs to.
%
% Optional input BINS is the number of bins used, and will instigate 
%   Gaussian smoothing of the tuning curves, fi(x).
%
%   In practice:
%       'sample' and 'training' are matrices with time samples for rows and
%           cell firing rates for columns.
%       'group' is a vector of spatial bins IDs describing where the mouse 
%            was at each 
%           corresponding time sample (row) of the matrix 'training'.
%
% CLASS
%  indicates which group each row of SAMPLE has been assigned to, and is
%  of the same type as GROUP.
%
% POSTERIOR, a matrix
%  containing estimates of the posterior probabilities that the j'th
%  training group was the source of the i'th sample observation, i.e.
%  Pr{group j | obs i}.
%
%   In practice:
%       class is a vector of absolute guesses (in the form of pixel ids) as
%           to where the rat was at each time sample (row) of the matrix 
%           'sample'.
%       posterior is a matrix of distributed guesses (in the form of a
%           probability distribution over all pixel ids) to where the rat 
%           was at each time sample (row) of the matrix 'sample'. Rows
%           correspond to time samples, columns correspond to pixel ids.
%
%

%% classify inputs
%
    if nargin > 4
        switch nargin
            case 5
                nbins_smoothing = varargin{1};
            otherwise
                error(message('too many inputs')) 
        end
    end

    
%check inputs
%
    if size(group,1) ~= size(training,1)
        error('stats:classify:InputSizeMismatch',...
            'The length of GROUP must equal the number of rows in TRAINING.');
    elseif size(sample,2) ~= size(training,2)
        error('stats:classify:InputSizeMismatch',...
            'SAMPLE and TRAINING must have the same number of columns.');
    end
 
%% Determine prior probabilities and Ca activity from the training data
%
% Here, we use the training data to calculate the prior probability of the 
%   occupying each state (spatial bin), and the activity rate maps of each cell over 
%   all states.

    %PRIOR PROBABILITY, P(S), i.e. the spatial occupancy distribution 
    %
    % P(S) = N(S) / sum(N(S))
    %
    % N(S) is the number of timesamples with the mouse in state (B)
    % 
    % This is the number of timesamples with the mouse in each bin 
    % divided by the total number of timesamples
    %  
    % This gives us a row vector length pixels
    %

% Group (bin) occupany
    [gindex,~,glevels] = grp2idx(group); % Returns a vector the length of the group at each step as well as a vector of the unique bins.
    gsize = histcounts(gindex); % Determine how many timesteps each bin was occupied
    
    occupancy_vector = (gsize/sum(gsize)); % P(S), spatial occupancy
    
    % Uniform prior
    occupancy_vector = ones(size(occupancy_vector))/length(occupancy_vector);
 
    %%   
 %FIRING RATE MAPS, f(x), i.e. tuning curve distribution
    %
    % f(x) is the mean firing rate observed in each pixel A for each cell
    %
    % This gives us a matrix bins*cells
    
    ngroups = numel(glevels); % Number of groups (spatial bins)
    ncells = size(training, 2); % The number of cells used
    
    % f_x = nan(size(ngroups, ncells)); % Preallocate a groups (spatial bins) x n cells matrix for cell activity
    
    
    % Check if traces are deconvolved spikes (s) or Calcium traces (C or
    % C_raw) and compute rate for spikes
    
    %% Calculate activity biases
    % Calculate P(A), the probability of a neuron being active AND
    % Calculate P(A|S), the probability of a neuron being active given the
    % behavioural state (spatial bin) i.e. each neuron's tuning curve
    
    tuning_curves = nan(ngroups, ncells); % Preallocate bins x neurons matrix
    if unique(training(training > 0)) == 1 % If binarized activity
        prob_active = sum(training)/length(training); % Number of frames with activity divided by total frames
        for igroup = 1:ngroups
            active_frames = sum(training(gindex == igroup, :)); % Calculate total active frames for each neuron in the current bin (igroup)
            frames_in_igroup = sum(gindex == igroup); % Calculate the total number of frames in the current bin
            tuning_curves(igroup,:) = active_frames./(frames_in_igroup);
        end        
    else % If flourescence compute mean
        prob_active = mean(training); % Mean flourescence activity over all frames
        for igroup = 1:ngroups
            %f_x(igroup,:) = mean(training(gindex == igroup,:)); % Calculate the mean flourescence intensity for a given group (spatial bin)
            tuning_curves(igroup,:) = mean(training(gindex == igroup,:)); % Calculate the mean flourescence intensity for a given group (spatial bin)
        end
    end
    
    %% Calculate the posterior
    % P(S|A) is the posterior probability distribution of behavioural state
    % (S) given the neural activity state (A). It is calculated as the
    % product of the observed likelihood of neural activity in an observed
    % behavioural state P(A|S) and the prior spatial distribution P(S)
    % divided by the bias of a neuron to fire P(A).
    
    % P_S and P_A are vectors corresponding to each cell, replicate them to
    % match the size of P_A_S, a bin x cell matrix
    occupancy_mtx = repmat(occupancy_vector', 1, ncells);
    prob_active_mtx = repmat(prob_active, ngroups, 1); 
       
    % P(A) can be smoothed with a Gaussian kernel
    %%%%%%%%% To be completed later %%%%%%%%%%%%%%%%%%%%%
    if exist('nbins_smoothing','var')
        kernel = gausswin(nbins_smoothing);
        kernel = kernel/sum(kernel);
        tuning_curves_smoothed = filter(kernel, 1, tuning_curves); % Filter applies the Gaussian kernel to each of the columns of f_x. I do not understand the second input
        tuning_curves = tuning_curves_smoothed;
    end
        


    %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%% FOR THE SAMPLE DATA %%%%%%%%%%%%%%%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%
% Here, we use the above calculated distributions to calculate the
%   conditional, P(B|A), and then the posterior, P(A|B), for the sample 
%   data.
%
    
    %CONDITIONAL PROBABILITY, P(B|A), i.e. the likelihood
    %
    % This is the probability for the firing rates B to occur given that 
    %   the animal is in pixel A. If we assume that the rates have 
    %  	Poisson distributions (are as random as possible in time) and 
    %  	that different cells are statistically independent of one another, 
    %  	then we can obtain the explicit expression:
    %
    %           P(B|A) = multiply_over P(Bi|A)
    %               and
    %           P(Bi|A) = multiply_over (((tao*fi(A))^Bi)/Bi!)*e^(-tao*fi(A))
    %
    % Here, P(B) is not explicitly calculated, but is obtained by
    % 	assuming that P(A|B) must equal 1.
    %
    % This gives us a matrix sample_timesamples*pixels
    %
    
   %% Decode the sample
   
   % For each timestep in calculate the conditional probability 
  
   n_sample_frames = size(sample,1); % The number of sample timesteps to classify is the number of rows
   
   decoded_probabilities = nan(ngroups, n_sample_frames); % Preallocate a matrix to hold the posterior probabilites for each testing sample timestep 
   bayesian_step_prob = nan(ngroups, ncells); % Preallocate a matrix to hold data during each loop
   
   for itimestep = 1:n_sample_frames
       active_cells = sample(itimestep,:) == 1; % Identify cells currently active in this frame
       inactive_cells = sample(itimestep,:) == 0; % Same for inactive cells
       bayesian_step_prob(:,active_cells) = (tuning_curves(:,active_cells).*occupancy_mtx(:,active_cells))./prob_active_mtx(:,active_cells);
       bayesian_step_prob(:,inactive_cells) = ((1 - tuning_curves(:,inactive_cells)).*occupancy_mtx(:,inactive_cells))./(1 - prob_active_mtx(:,inactive_cells));
       
       decoded_probabilities(:,itimestep) = expm1(sum(log1p(bayesian_step_prob),2)); % This should be used instead of simple product to avoid numerical underflow
       decoded_probabilities(:,itimestep) = decoded_probabilities(:,itimestep)./sum(decoded_probabilities(:,itimestep),'omitnan'); % The method above is an approximation, to obtain a posterior distribution, we can normalize every datapoint
    end
       
   %% Outputs
   posterior = decoded_probabilities;
   [~, class] = max(decoded_probabilities);
   class = class'; % Output as a column vector
    











