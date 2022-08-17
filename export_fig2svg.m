function export_fig2svg(path)
% Function that takes a folderpath and converts the matlab .fig files to
% .svg files

if exist(path) == 7 % check if received a folder
    folderpath = path;

    contents = dir(folderpath); % Get a listing of the files in the folder
    contents = contents(contains({contents.name},'.fig')); % Drop any files found that aren't .fig
    
    if isempty(contents)
        error('No .fig files found in folderpath')
    end
    

    for i = 1:numel(contents)
        fig_file = fullfile(contents(i).folder, contents(i).name); % Get the full file path of each .fig file
        current_figure = openfig(fig_file); % open the figure
        svg_file = strrep(fig_file,'.fig','.svg'); % Create the output .svg filepath in the same folder
        print(current_figure, svg_file,'-dsvg','-painters') % Use print to export. dsvg specifies file type, painters specifices editable vector file
        close(current_figure);
        disp(['Exported ', contents(i).name, ' to .svg file'])
    end
elseif exist(path) == 2 % check if received a file
    fig_file = path;
    if ~contains(fig_file,'.fig')
        error('File is not a .fig file.')
    end
    svg_file = strrep(fig_file,'.fig','.svg'); % Create the output .svg filepath in the same folder
    current_figure = openfig(fig_file); % open the figure
    print(current_figure, svg_file,'-dsvg','-painters') % Use print to export. dsvg specifies file type, painters specifices editable vector file
    close(current_figure);
    [~, file_name, file_ext] = fileparts(fig_file);
    disp(['Exported ', fullfile(file_name,file_ext), ' to .svg file'])
else
    error('Input path is not a folder nor a filepath.')
end


