function [behavior_mtx, traces] = load_mitchdata(foldername)
%loads data from a CNMFE and deepLab cut output containing folder

%behavior_mtx key
% col#  data
%   1   time
%   2   xpos
%   3   ypos
%   4   trial number
%   5   velocity
%   6   acceleration



% TO DO
% plot delete-jumps function effects to determine best parameters



%% Load CNMFE output (traces and time)
[traces, frame_times] = load_cnmfe(foldername);



%% Load DEEP LAB CUT output (camera position)

% Load position time
dlc_time_file = get_file_paths_targeted(foldername, {'dlc_time','.csv'});
pos_times = readmatrix(dlc_time_file{1});

% load position deep lab cut xy matrix (has two different naming
% conventions, dending on computer used)
dlc_xy_file = get_file_paths_targeted(foldername, {'behaviourDeepCut','.csv'});
dlc_xy_file = [dlc_xy_file; get_file_paths_targeted(foldername, {'behaviourDLC','.csv'})];
dlc_out = readmatrix(dlc_xy_file{1});

% xy head positions
head_xy = dlc_out(:, [5 6]); % head position (x,y)
head_xy(dlc_out(:,7)<0.99, :) = nan; % nan low confidence positions

% combine time and position
behavior_mtx = [pos_times head_xy];

% cup positions
leftCup_xy = nanmean(dlc_out(dlc_out(:,13)>0.99, [11 12])); % left cup
rightCup_xy = nanmean(dlc_out(dlc_out(:,22)>0.99, [20 21])); % right cup



%% Movie frame overlay
%{
figure; hold on

try
    mov_file = get_file_paths_targeted(foldername, {'.avi'});
    mov = VideoReader(mov_file{1});
    im = readFrame(mov);
    image(im);
catch
end

plot(head_xy(:,1), head_xy(:,2), 'y-')
plot(leftCup_xy(1), leftCup_xy(2), '.', 'markersize', 50)
plot(rightCup_xy(1), rightCup_xy(2), '.', 'markersize', 50)
%}



%% Resample time

% catch computer glitch where different number of calcium frames than frame
% times
if size(traces,2) > size(frame_times,1)
    warning(['glitch; fewer frame times than calcium frames for ', foldername])
    traces = traces(:,1:size(frame_times,1),:);
elseif size(traces,2) < size(frame_times,1)
    warning(['glitch; more frame times than calcium frames for ', foldername])
    frame_times = frame_times(1:size(traces,2),:,:);
end

% find common time
min_max_flor_time = [min(frame_times) max(frame_times)];
min_max_beh_time = [min(behavior_mtx(:,1)) max(behavior_mtx(:,1))];
min_max_common_time = [max([min_max_flor_time(1) min_max_beh_time(1)]) min([min_max_flor_time(2) min_max_beh_time(2)])];

% cut edges of session to align time
behavior_mtx(behavior_mtx(:,1)<min_max_common_time(1) | behavior_mtx(:,1)>min_max_common_time(2), :) = [];
traces(:, frame_times<min_max_common_time(1) | frame_times>min_max_common_time(2), :) = [];
frame_times(frame_times<min_max_common_time(1) | frame_times>min_max_common_time(2)) = [];

% resample behavior at 100hz
frame_times = frame_times - min(behavior_mtx(:,1)); % Zero the Ca trace timestamps
behavior_mtx(:,1) = behavior_mtx(:,1) - min(behavior_mtx(:,1)); % Zero the behaviour timestamps
new_frame_rate = 20; % 20 Hz
behavior_mtx = resample_time(behavior_mtx,new_frame_rate); % MDS changed from 100 Hz to 20 Hz 2022-05-13


% preallocate resampled flor at 100hz
new_flor = zeros(size(traces,1), size(behavior_mtx,1), 3);

% upsample flor
for icell = 1:size(traces,1)
    
    % load resampled S traces 
    % place spike in new 100hz time bin that is closest to old 20hz time bin
    all_spike_times = frame_times(traces(icell,:,1)==1);
    new_times_temp = repmat(behavior_mtx(:,1),[1 length(all_spike_times)]);
    [~,closestIndex] = min(abs(new_times_temp-all_spike_times'));
    new_flor(icell,closestIndex,1) = 1; 
    
    % load resampled C and C_raw traces   
    % First and last few timestamp columns may become NaNs because the
    % frame rates of behaviour and calcium cameras are different so zeroing
    % by the behaviour camera doesn't necessarily but the trace time to
    % zero. Could maybe use 'extrap' marker to incorporate those times.
    new_flor(icell,:,2) = interp1(frame_times, traces(icell,:,2), behavior_mtx(:,1), 'linear');
    new_flor(icell,:,3) = interp1(frame_times, traces(icell,:,3), behavior_mtx(:,1), 'linear');
    
end
traces = new_flor; clearvars new_flor;



%% Correct position

% rotate positions
[behavior_mtx] = rotate_pos(leftCup_xy, rightCup_xy, behavior_mtx);
[pos_out] = rotate_pos(leftCup_xy, rightCup_xy, [nan leftCup_xy; nan rightCup_xy]);
leftCup_xy = pos_out(1,2:3); rightCup_xy = pos_out(2,2:3);


%figure; hold on
%patch([leftCup_xy(1)-20 rightCup_xy(1)+20 rightCup_xy(1)+20 leftCup_xy(1)-20], [leftCup_xy(2)-30 leftCup_xy(2)-30 leftCup_xy(2)+20 leftCup_xy(2)+20],'r')
%plot(behavior_mtx(:,2), behavior_mtx(:,3), 'k.'); title hard-cut

%timeidx = behavior_mtx(:,1)>350 & behavior_mtx(:,1)<450;
%figure; plot3(behavior_mtx(timeidx,2), behavior_mtx(timeidx,3), behavior_mtx(timeidx,1)); title orig



% delete out of range items % MDS can come in here and use the measurements
% from the sensors to produce a more precise rectangle
accepted_area = [leftCup_xy(1)-20 rightCup_xy(1)+20 leftCup_xy(2)-30 leftCup_xy(2)+13]; % 13 has tuned down from 20 to account for streaks.
behavior_mtx = hard_correct_pos(behavior_mtx, accepted_area(1:2), accepted_area(3:4));


%figure; plot3(behavior_mtx(:,1), behavior_mtx(:,2), behavior_mtx(:,3), 'k-'); title post-hard-cut

% remove implausible position changes
behavior_mtx = resample_jumps(behavior_mtx, 800, 200);

%figure; plot3(behavior_mtx(timeidx,2), behavior_mtx(timeidx,3), behavior_mtx(timeidx,1)); title delete

% interpolate nans
behavior_mtx = interp_at_nans(behavior_mtx);

%figure; plot3(behavior_mtx(timeidx,2), behavior_mtx(timeidx,3), behavior_mtx(timeidx,1)); title interp

% smooth positions
behavior_mtx = smooth_pos(behavior_mtx, 9);


%figure; plot3(behavior_mtx(timeidx,2), behavior_mtx(timeidx,3), behavior_mtx(timeidx,1)); title smooth


%% Normalize position using reward locations

% subtract
behavior_mtx(:,2) = behavior_mtx(:,2)-leftCup_xy(1);
behavior_mtx(:,3) = behavior_mtx(:,3)-leftCup_xy(2);

% reward cup x position difference
cup_diff = (rightCup_xy(1)- leftCup_xy(1));

%normalize
behavior_mtx(:,2) = behavior_mtx(:,2) ./ cup_diff;
behavior_mtx(:,3) = behavior_mtx(:,3) ./ cup_diff;

% update reward locations
rightCup_xy(1) = (rightCup_xy(1)-leftCup_xy(1))./cup_diff;
leftCup_xy(1) = leftCup_xy(1) - leftCup_xy(1);
rightCup_xy(2) = rightCup_xy(2) - leftCup_xy(2);
leftCup_xy(2) = leftCup_xy(2) - leftCup_xy(2);

%figure; image(im); hold on
%
figure; hold on
plot(behavior_mtx(:,2), behavior_mtx(:,3), 'r-')
plot(leftCup_xy(1), leftCup_xy(2), '.', 'markersize', 50)
plot(rightCup_xy(1), rightCup_xy(2), '.', 'markersize', 50)
maze_outline
%}



%% Trials

%use reward areas as ends of maze
min_max_x = [leftCup_xy(1) rightCup_xy(1)];
x_end_bounds = linspace(min_max_x(1), min_max_x(2), 21); % 2022-05-09 MDS changed 10 bins to 20 to try only cutting out 5% each end
end_lo = x_end_bounds(2);
end_hi = x_end_bounds(end-1);

%iterate through time points identifying trials
trial_label = nan(size(behavior_mtx,1),1);
previous_sect = 0;
cross_times = [];
trial_count = 0;

%find trial boundary crossings
for itime = 1:size(behavior_mtx,1)

    if isnan(behavior_mtx(itime,2))
        continue
    end
    
    %if coming from center
    if ismember(previous_sect, [0 2])       
        current_sect = find_current_sect(behavior_mtx(itime,2), end_lo, end_hi);
    %if coming from end   
    elseif ismember(previous_sect, [0 1 3])
        current_sect = find_current_sect(behavior_mtx(itime,2), end_lo, end_hi);
    end
    
    %crossing
    if previous_sect~=0 && previous_sect~=current_sect
        sections = [previous_sect current_sect];
        current_end = sections(ismember(sections, [1 3]));
        cross_times = [cross_times; [behavior_mtx(itime,1) current_end]];
    end
    
    %shift forward
    previous_sect=current_sect;
end

%number trials 
while 1
    trial_count = trial_count+1;
    first1 = min(cross_times(cross_times(:,2)==1,1));
    first3 = min(cross_times(cross_times(:,2)==3,1));
    trial_end = max([first1 first3]);
    if first1 < first3
        trial_start = max(cross_times(cross_times(:,2)==1 & cross_times(:,1)<trial_end,1));
        %start_side = 1;
    elseif first1 > first3
        trial_start = max(cross_times(cross_times(:,2)==3 & cross_times(:,1)<trial_end,1));
        %start_side = 3;
    else
        break
    end
    
    %load
    trial_label(behavior_mtx(:,1)>=trial_start & behavior_mtx(:,1)<trial_end) = trial_count;
    
    %delete used crossings
    cross_times = cross_times(cross_times(:,1)>trial_end, :);
    
end

% add trial labels to behavior matrix
behavior_mtx = [behavior_mtx trial_label];

% trial time bounds
max_trialduration = 10; %s
num_trials = length(unique(trial_label(~isnan(trial_label))));
trial_time_bounds = nan(num_trials,2);
for itrial = 1:num_trials 
    trial_time_bounds(itrial,:) = [nanmin(behavior_mtx(trial_label==itrial,1)) nanmax(behavior_mtx(trial_label==itrial,1))];
end
trial_durations = diff(trial_time_bounds');

% exclude trials
unq_trials = unique(behavior_mtx(~isnan(behavior_mtx(:,4)),4));
excluded_trials = unq_trials(trial_durations > max_trialduration);

% plot trial paths and durations
figure; set(gcf, 'Position', [25 25 1065 476]) % MDS lowered from 677 373
mouse_session_str_parts = regexp(foldername,filesep,'split');
%mouse_session_string = sprintf(
sgtitle(['Trial behavior ',mouse_session_str_parts{end-1},' ',mouse_session_str_parts{end}]);
subplot(1,2,1); hold on; 
plot3(behavior_mtx(:,2), behavior_mtx(:,1), behavior_mtx(:,3), '-', 'color', 0.7.*[1 1 1])
for itrl = unique(trial_label)'
    trl_x = behavior_mtx(trial_label==itrl,2);
    trl_t = behavior_mtx(trial_label==itrl,1);
    trl_y = behavior_mtx(trial_label==itrl,3);
    
    % excluded trials plotted with black
    if ismember(itrl,excluded_trials)
        plot3(trl_x, trl_t, trl_y, 'k-', 'linewidth', 2)
    else
        plot3(trl_x, trl_t, trl_y, 'linewidth', 2)
    end
en
ylim([-10 inf])
ylabel('time (s)')
xlabel('x position')
set(gca,'TickLength',[0, 0]); box off;

subplot(1,2,2); hold on
histogram(trial_durations, 0:10:200); 
plot(max_trialduration.*[1 1], ylim, 'r-')
ylabel('trial count')
xlabel('duration (s)')
set(gca,'TickLength',[0, 0]); box off;

% update behavior matrix
behavior_mtx(ismember(behavior_mtx(:,4), excluded_trials),4) = nan;
for itrl = 1:length(unq_trials)
    current_trial = unq_trials(itrl);
    behavior_mtx(behavior_mtx(:,4)==current_trial,4) = itrl;
end



%% Velocity

velocity_vect = nan(size(behavior_mtx,1),1);
time_jump = 21;
for i = 1:length(behavior_mtx(:,2))-time_jump
    %pos points for evaluating velocity
    p1 = behavior_mtx(i, 2:3);
    p2 = behavior_mtx(i+time_jump, 2:3);

    %time points for evaluating velocity
    t1 = behavior_mtx(i,1);
    t2 = behavior_mtx(i+time_jump,1);
    
    %velocity
    velocity_vect(i+10) = pdist([p1; p2])/(t2-t1);
end
maze_length_in_pixels = pdist([leftCup_xy;rightCup_xy]);
maze_length_in_cm = 91;
pixels_per_cm = maze_length_in_pixels/maze_length_in_cm; % Error is here, maze position has been normalized to 1 already
pixels_per_cm_over_time_jump = pixels_per_cm*time_jump;
behavior_mtx(:,5) = velocity_vect;
behavior_mtx(:,5) = behavior_mtx(:,5).*pixels_per_cm_over_time_jump; %cm/s
%behavior_mtx(:,5) = behavior_mtx(:,5)./100; %m/s



%% Acceleration

acceleration_vect = nan(size(behavior_mtx,1),1);
for i = 1:length(behavior_mtx(:,2))-(time_jump+1)
    
    %vel points for evaluating velocity
    p1 = velocity_vect(i);
    p2 = velocity_vect(i+time_jump);

    %time points for evaluating velocity
    t1 = behavior_mtx(i,1);
    t2 = behavior_mtx(i+time_jump,1);
    
    %velocity
    acceleration_vect(i+10) = (p2-p1)/(t2-t1);
end
behavior_mtx(:,6) = acceleration_vect;



end




%% Internal functions

function timeXY = hard_correct_pos(timeXY, lohiX, lohiY)
%nan any pos values outside of rectangle defined by lohiX and lohiY
outside_rectX_idx = timeXY(:,2)<lohiX(1) | timeXY(:,2)>lohiX(2); 
outside_rectY_idx = timeXY(:,3)<lohiY(1) | timeXY(:,3)>lohiY(2);

%nan points
timeXY(outside_rectX_idx | outside_rectY_idx, [2 3]) = nan;

end

function [timeXY, dists] = resample_jumps(timeXY, too_big_origin, cum_mod)
%delete eroneous spatial positions

    %prep to remove improbable changes in position
    too_big = too_big_origin;
    guilty_by_association = 4;

    %adjacent pos points for evaluating velocity
    first_nnan_row = find(~isnan(timeXY(:,2)),1,'first');
    p1 = timeXY(first_nnan_row, 2:3);
    p2 = timeXY(first_nnan_row+1, 2:3);

    %adjacent time points for evaluating velocity
    t1 = timeXY(first_nnan_row,1);
    t2 = timeXY(first_nnan_row+1,1);

    %preallocate
    deletions = zeros(size(timeXY,1),1);
    dists = zeros(size(timeXY,1),1);

    %iterate through adjacent points to evaluate velocity
    count = 0;
    for i = first_nnan_row:length(timeXY(:,2))-2

        %velocity
        current_distance = pdist([p1; p2])/(t2-t1);
        dists(i+1) = current_distance;

        %if the current velocity is too big
        if ~(current_distance <= too_big)

            %note that point (and the next 4) should be deleted (index for later)
            if length(timeXY(:,2))-2-i > guilty_by_association
                deletions(i:i+guilty_by_association) = 1;
            end

            %move to the next point, but keep the first of the adjacent pair
            p2 = timeXY(i+2, 2:3);
            t2 = timeXY(i+2, 1);

            %each time it's too big, increase what is considered "too big"
            count = count + cum_mod;
            too_big = too_big + count;

        %if it's not too big
        else

            %reset what is considered "too big"
            too_big = too_big_origin;
            count = 0;

            %update points
            p1 = timeXY(i+1, 2:3);
            p2 = timeXY(i+2, 2:3);
            t1 = timeXY(i+1, 1);
            t2 = timeXY(i+2, 1);

        end
    end

    
    %figure; histogram(dists,0:100:3000); ylim([0 100])
    
    %index to delete dubious points
    deletions = logical(deletions);
    
    %{
    timeidx = timeXY(:,1)>350 & timeXY(:,1)<450;
    figure; plot3(timeXY(timeidx,2), timeXY(timeidx,3), timeXY(timeidx,1), 'k-'); title orig
        hold on; plot3(timeXY(timeidx & ~deletions,2), timeXY(timeidx & ~deletions,3), timeXY(timeidx & ~deletions,1), 'b.');
    hold on; plot3(timeXY(timeidx & deletions,2), timeXY(timeidx & deletions,3), timeXY(timeidx & deletions,1), 'r.');
    %}
    
    timeXY(deletions, 2:3) = NaN;

    
    
    %figure; plot3(timeXY(timeidx,2), timeXY(timeidx,3), timeXY(timeidx,1)); title after

end

function timeXY = resample_time(timeXY,new_frame_rate)
    %resample time and interpolate new position values
    mtx(:,1) = (0:(1/new_frame_rate):max(timeXY(:,1)))';
    mtx(:,2) = interp1(timeXY(:,1), timeXY(:,2), mtx(:,1), 'linear');
    mtx(:,3) = interp1(timeXY(:,1), timeXY(:,3), mtx(:,1), 'linear');
    timeXY = mtx;
end

function timeXY = interp_at_nans(timeXY)
%replace nan postions with interpolated values
    non_nan_pos = ~isnan(timeXY(:,2)) & ~isnan(timeXY(:,3)); %pos index
    timeXY(:,2) = interp1(timeXY(non_nan_pos, 1), timeXY(non_nan_pos, 2), timeXY(:,1), 'linear');
    timeXY(:,3) = interp1(timeXY(non_nan_pos, 1), timeXY(non_nan_pos, 3), timeXY(:,1), 'linear');
end

function timeXY = smooth_pos(timeXY, smoothing_window_size)
    non_nan_pos = ~isnan(timeXY(:,2)) & ~isnan(timeXY(:,3)); %index
    timeXY(non_nan_pos,2) = smooth(timeXY(non_nan_pos, 2), smoothing_window_size);
    timeXY(non_nan_pos,3) = smooth(timeXY(non_nan_pos, 3), smoothing_window_size);
end

function [timeXY, rotang] = rotate_pos(leftAnchor_xy, rightAnchor_xy, timeXY)
%rotates positions and and HDs so the track (mouse's path) is parallel with
%the x axis

    % compute counterclockwise rotation angle from two anchor points
    rotang = -atan2(rightAnchor_xy(2) - leftAnchor_xy(2), rightAnchor_xy(1) - leftAnchor_xy(1));
    
    %build rotation matrix
    romat=[cos(rotang) -sin(rotang);sin(rotang) cos(rotang)];       

    %rotation occurs around origin, so we temporarily center all points around the origin
    timeXY(:,2) = timeXY(:,2) - mean([leftAnchor_xy(1) rightAnchor_xy(1)]);
    timeXY(:,3) = timeXY(:,3) - mean([leftAnchor_xy(2) rightAnchor_xy(2)]);

    %apply rotation
    timeXY(:,[2 3]) = (romat*timeXY(:,[2 3])')';

    %undo centering
    timeXY(:,2) = timeXY(:,2) + mean([leftAnchor_xy(1) rightAnchor_xy(1)]);
    timeXY(:,3) = timeXY(:,3) + mean([leftAnchor_xy(2) rightAnchor_xy(2)]);      
    
end

function current_sect = find_current_sect(trial_pos_x, end_lo, end_hi)
    %if on the center
    if trial_pos_x>end_lo && trial_pos_x<end_hi
        current_sect = 2;
    %if on L
    elseif trial_pos_x<end_lo
        current_sect = 1;
    %if on R
    elseif trial_pos_x>end_hi
        current_sect = 3;
    end
end
