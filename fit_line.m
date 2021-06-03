function [R, p, se, poly] = fit_line(X,Y, varargin)
%fits a line to a scatterplot of X,Y
% plot style
% 0 = line only
% 1 = dots only
% 2 = line and dots

% input
if nargin == 3
    plot_style = varargin{1};
else
    plot_style = 2;
end

poly=polyfit(X(~isnan(X) & ~isnan(Y)),Y(~isnan(X) & ~isnan(Y)),1);
fit_x = (min(X)-(range(X)/4)):(range(X)/100):(max(X)+(range(X)/4));
fit_y=polyval(poly, fit_x);

% square error
se = (Y-polyval(poly, X)).^2;

%figure;
hold on; 
if plot_style == 0
    plot(fit_x, fit_y, 'k-', 'LineWidth', 2)
elseif plot_style == 1
    plot(X,Y,'ko');
elseif plot_style == 2
    plot(X,Y,'ko');
    plot(fit_x, fit_y, 'k-', 'LineWidth', 2)
end


[R, p] = corr(X(~isnan(X) & ~isnan(Y)),Y(~isnan(X) & ~isnan(Y)));

set(gca,'TickLength',[0, 0]); box off; %axis square

end