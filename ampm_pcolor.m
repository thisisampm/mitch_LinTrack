function ampm_pcolor(mtx_in)
% uses pcolor in the fashion of imagesc

% pad sides
mtx_in = [mtx_in zeros(size(mtx_in,1),1) ];
mtx_in = [mtx_in; zeros(1, size(mtx_in,2))];

% plot
h = pcolor(mtx_in);

% no black lines
set(h, 'EdgeColor', 'none');
set(gca,'TickLength',[0, 0]);

% set yaxis
set(gca, 'YDir','reverse')

% center tick marks
%{
xticks(1.5:1:size(mtx_in,2))
xticklabels(xticks-0.5)
yticks(1.5:1:size(mtx_in,1))
yticklabels(yticks-0.5)
%}

