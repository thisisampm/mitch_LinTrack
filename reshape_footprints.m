%% Reshape cnfme spatial outputs for use in CellReg
% Take a mouse directory with subdirectories containing CNMFE outputs organized by day
% and reshapes the footprints for CellReg

function [] = reshape_footprints(mouse_folder)

%% Get all the cnmfe files 

% mouse_folder_full_path = fullfile('C:\Users\mitch\OneDrive - University of Toronto\PhD\Experiments\LinearTrack',mouse_folder_number);
files = dir(fullfile(mouse_folder,'*','*_cnmfe.mat'));

%% If necessary, create folder to store footprints.
savedir = fullfile(mouse_folder,'Footprints');

if ~exist(savedir)
    mkdir(savedir);
    disp(['Footprints folder created: ', savedir]);
end
   
%% Loop through the cnmfe files reshaping the spatial footprints
for i = 1:numel(files)
    % Load the current Ca file
    Cafile = [files(i).folder, '/', files(i).name];
    disp(['Reshaping: ', files(i).name])
    warning('off'); load(Cafile); warning('on');
    % Preallocate a 3d matrix to store footprints neurons x y pixels x
    % xpixels
    fp = nan([size(st.A,2), size(st.Cn)]);
    % For each neuron, change the spatial footprint from vector to matrix
    for j = 1:size(st.A,2)
       fp(j,:,:) =  reshape(st.A(:,j),size(st.Cn));
    end
    % Create a pathname for each session's footprints from output
    savename = [savedir,'/', erase(files(i).name, '_cnmfe.mat'), '_footprints.mat'];
    % Save the footprints
    save(savename, 'fp')
    %clear st fp
end