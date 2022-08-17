
function errorbar_scatter(vect_cell, varargin)
% overlays jittered plots of each vector in vect_cell to means xpositions
% 1 : length(vect_cell)

% connect optional input
if nargin == 2
    connect = varargin{1};
    xaxis_input = 1:length(vect_cell);
    colors = linspecer(length(vect_cell));
elseif nargin == 3
    connect = varargin{1};
    xaxis_input = varargin{2};
    colors = linspecer(length(vect_cell));
elseif nargin == 4
    connect = varargin{1};
    xaxis_input = varargin{2};
    colors = varargin{3};
else
    connect = 0;
    xaxis_input = 1:length(vect_cell);
    colors = linspecer(length(vect_cell));
end


% grey or color
max_vect = nan(size(vect_cell)); 
for ivect = 1:length(vect_cell); max_vect(ivect) = length(vect_cell{ivect}); end 
max_vect = nanmax(max_vect);
%colors = repmat(0.8.*[1 1 1], length(vect_cell),1); % gray
%colors = distinguishable_colors(max_vect);

hold on
jxpos = cell(1,length(vect_cell));
for ivect = 1:length(vect_cell)
    jxpos{ivect} = jitter_xpos(xaxis_input(ivect), vect_cell{ivect});
    for id = 1:length(vect_cell{ivect})
        plot(jxpos{ivect}(id), vect_cell{ivect}(id),'o', 'color', colors(id,:))
    end
end

if connect == 1
    for ivect = 1:length(vect_cell)
        vect_cell{ivect} = vect_cell{ivect}(:);
        jxpos{ivect} = jxpos{ivect}(:);
    end
    vect_mtx = cell2mat(vect_cell);
    jxpos = cell2mat(jxpos);
    
    for ivect = 1:size(vect_mtx,1)
        plot(jxpos(ivect,:), vect_mtx(ivect,:), '-', 'color', colors(ivect,:))
    end
end

