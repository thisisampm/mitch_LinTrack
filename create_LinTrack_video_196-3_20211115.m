%%
%% Load data
% Load the calcium video, behaviour video and cnmfe with traces and
% contours

folder_path = 'C:\Users\mitch\OneDrive - University of Toronto\PhD\Experiments\LinearTrack\Raw_track_data\control\196-3';

Ca_vid_path = fullfile(folder_path, '2021-11-25-10-16-35.mkv');
cnmfe_path = fullfile(folder_path, '196-3_2021-11-25-10-16-35_cnmfe.mat');
pos_vid_path = fullfile(folder_path, '2021_11_25_10_16_33_behaviour.avi');


if ~exist('st')
    load(cnmfe_path)
end

% Load matching timestamps
behav_ts_path = fullfile(folder_path, '2021_11_25_10_16_33_ts_dlc_time.csv');
Ca_ts_path = fullfile(folder_path, '2021-11-25-10-16-35_Ca_time.csv');

behav_ts = table2array(readtable(behav_ts_path));
Ca_ts = table2array(readtable(Ca_ts_path));

% Access video files
Ca_v = VideoReader(Ca_vid_path);
pos_v = VideoReader(pos_vid_path);

%% Read desired video frames
% For this video want 0:51-1:18 from behaviour video timestamps
start_time = 45; 
end_time = 78; 
pos_start = floor(start_time*pos_v.FrameRate);
pos_end = ceil(end_time*pos_v.FrameRate);
pos_frames = read(pos_v,[pos_start pos_end]);
ROI = (190:290); % Define an roi to crop only around the track
pos_resized = imresize(pos_frames(ROI,:,:,:),(Ca_v.Width*2/pos_v.Width)); % Rescale the behavioural video to be twice the width of the Ca video


% Ca video is at a different frame rate than behaviour video, therefore get
% closest frames to behaviour video timestamps and read the frames in
% between
[~, ca_start] = min(abs(Ca_ts - behav_ts(pos_start)));
[~,ca_end] = min(abs(Ca_ts - behav_ts(pos_end)));
ca_frames = read(Ca_v, [ca_start ca_end]);

[~, closest_Ca] = min(abs(behav_ts(pos_start:pos_end) - Ca_ts(ca_start:ca_end)'),[],2);


%% Preallocate structure to hold output video frames
nFrames = size(pos_resized,4); % Number of frames equal to the number of behavioural video frames displayed
s(nFrames) = struct('cdata',[],'colormap',[]); % Structure to hold movie data from getframe function. Must be in this format

%% Set up figure, axes to hold Ca image, axes to hold behaviour image, and plot Ca traces
% This allows us to keep everything positioned in the proper order, and
% just change the image data and the plot location every loop iteration
hFig = figure('MenuBar','none',... % Whole figure
    'Units','pixels',...
    'Position',[4000 -300 (Ca_v.Width*2) (Ca_v.Height + size(pos_resized,1) + 10)]);
hAx_coor = axes('Parent',hFig,... % Plot the contours and hold Ca video image data
    'Units','pixels',...
    'Position',[0 0 Ca_v.Width Ca_v.Height],...
    'NextPlot','add',...
    'Visible','off',...
    'XTick',[],...
    'YTick',[]);
hIm_Ca = imshow(uint8(zeros(Ca_v.Height,Ca_v.Width,3)),... % Ca video image data
    'Parent',hAx_coor);
hAx_trace = axes('Parent',hFig,... % Axis to plot calcium traces
    'Units','pixels',...
    'Position',[Ca_v.Width 0 Ca_v.Width Ca_v.Height],...
    'NextPlot','add',...
    'Visible','off',...
    'XTick',[],...
    'YTick',[]);
hAx_pos = axes('Parent',hFig,... % Axis to hold the behavioural position video frame
    'Units','pixels',...
    'Position',[0 (Ca_v.Height+5) size(pos_resized,2) size(pos_resized,1)],...
    'NextPlot','add',...
    'Visible','off',...
    'XTick',[],...
    'YTick',[]);
hIm_pos = imshow(uint8(zeros(size(pos_resized,1), size(pos_resized,2),3)),... % Behaviour position video image data
    'Parent',hAx_pos);
axis(hAx_coor,'image');
axis(hAx_pos,'image'); % Change axes with images so that the image fills them completely

%% Get relevant cnmfe data and plot the contours and trace data
% Define and plot cell contours
%candidate_cells = [24,110,16,21,90,7]; % candidate_cells_196_3 = [181, 167, 145, 136, 134, 123, 116, 115, 111, 105, 104, 100, 95, 93, 89, 84, 80, 70, 69, 65, 64, 62, 60, 55, 54, 51, 50, 48, 46, 43, 42, 39, 32, 26, 22, 18, 17, 14, 12, 10, 7, 5, 4, 1]
%candidate_cells = [12, 3, 103, 107, 110, 118, 120, 124, 4, 127, 8, 136, 7, 9, 79, 82, 75, 60 ,61, 35,29];

candidate_cells = [1 ,3, 139, 18, 32, 30, 67, 45, 11, 26, 6, 8];
%candidate_cells = pc_ind;

coor = st.Coor(candidate_cells); % Get the contours of candidate cells
coor_color = distinguishable_colors(numel(coor)); % Define a unique color for each contour
% Plot each contour 
for i = 1:numel(coor)
    plot(hAx_coor,coor{i}(1,1:end),coor{i}(2,1:end),'color',coor_color(i,:),'LineWidth',1.25)
end

% Get Ca traces over the selected frames
C = st.C(candidate_cells, ca_start:ca_end); % Deconvoled traces
C = C - min(C,[],2); % Scale the traces
C = C./max(C,[],2);
C = C';
offset = linspace(0,size(C,2)/3,size(C,2)); % Create a vector to add to the traces causing them to be offset by a 3rd when plotting
C_offset = C + offset;

% Define the axes for plotting traces
xmax = size(C_offset, 1); % x axis is the length of time x cell traces matrix
ymax = max(max(C_offset)); % y axis height is the greatest trace height when considering the offset
hAx_trace.XLimMode = 'manual'; hAx_trace.YLimMode = 'manual';
axis(hAx_trace,[(0 - 0.05*xmax) (xmax + 0.05*xmax) (0 - 0.05*ymax) (ymax + 0.05*ymax) ]); % Set the axes with %5 padding around maximum and min values
hold(hAx_trace, 'on');  

% Plot the Ca traces
hLines = plot(hAx_trace, C_offset); % Plot each trace on the bottom right trace axis
coor_color_cell = mat2cell(coor_color, ones(1,size(coor_color,1))); % Break the color matrix into a cell array so that it works with plot
set(hLines,{'Color'},coor_color_cell,'LineWidth', 0.75); % Set line color using cell of colors

%% Loop through creating 

for i = 1:size(pos_resized,4)
    hIm_pos.CData = pos_resized(:,:,:,i); % Update behavioural video frame
    hIm_Ca.CData = ca_frames(:,:,:,closest_Ca(i)); % Update Ca video frame
    
    xdata = 1:closest_Ca(i);
    ydata = mat2cell(C_offset(xdata,:)',ones(1,size(C_offset,2))); % Create a cell array for each line with the data up to this frame
    set(hLines,{'Xdata'},{xdata},{'Ydata'},ydata); % Update the traces data being plotted
    %axis(hAx_trace,[0 xmax 0 ymax]);
    drawnow
        
    s(i) = getframe(hFig); % Capture each frame of entire figure.
end
hold(hAx_trace,'off');

%% Write movie to file

output_video_path = fullfile(folder_path, '2022-03-15_LinTrack_Example.mp4');

vout = VideoWriter(output_video_path, 'MPEG-4');
vout.FrameRate = 45; % Set framerate to 2x speed
open(vout);
writeVideo(vout,s) % Write the movie structure to the file
close(vout);
