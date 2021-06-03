function trace_footprint(matrix, varargin)
% plots outline of a single cell's spatial footprint
% input is a 2d matrix

% input
if length(varargin)==1
    plot_color = varargin{1};
    xpos = 0;
    ypos = 0;
    line_width = 1;
elseif length(varargin)==2
    plot_color = varargin{1};
    xpos = varargin{2};
    ypos = 0;
    line_width = 1;
elseif length(varargin)==3
    plot_color = varargin{1};
    xpos = varargin{2};
    ypos = varargin{3};
    line_width = 1;
elseif length(varargin) == 4
    plot_color = varargin{1};
    xpos = varargin{2};
    ypos = varargin{3};
    line_width = varargin{4};
else
    plot_color = [];
    xpos = 0;
    ypos = 0;
    line_width = 1;
end


% smooth footprint
matrix = double(matrix);
matrix = smooth2a(matrix, 5);
matrix = double(matrix>range(matrix(:))/2);

% aesthetics
hold on
ylim([0.5 size(matrix,1)+0.5])
xlim([0.5 size(matrix,2)+0.5])
set(gca, 'YDir','reverse')
set(gca,'TickLength',[0, 0]); box off;
xticks([])
yticks([])
axis off
axis square

%find contour coords using contour function
h = contourc(matrix,1);
h(:,1) = [];
end_field = find(h(1,:)==.5, 1, 'first');
h(:, end_field:end) = [];

%plot contour
%
smooth_wndw = 10;
if ~isempty(plot_color)
    plot(smooth([h(1,:) h(1,1)],smooth_wndw)+xpos, smooth([h(2,:), h(2,1)],smooth_wndw)+ypos, 'color', plot_color, 'Linewidth', line_width)
else
    plot(smooth([h(1,:) h(1,1)],smooth_wndw)+xpos, smooth([h(2,:), h(2,1)],smooth_wndw)+ypos, 'Linewidth', line_width)
end
%}
%plot([h(1,:) h(1,1)], [h(2,:), h(2,1)], 'color', plot_color)


end

            
