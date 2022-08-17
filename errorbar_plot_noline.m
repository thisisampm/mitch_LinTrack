function errorbar_plot_noline( cell_e , varargin)
%plots errorbar


%connections
if nargin == 2
    connect = varargin{1};
elseif nargin ==3
    connect = varargin{1};
    xaxis = varargin{2};
elseif nargin==4
    connect = varargin{1};
    xaxis = varargin{2};
    colors = repmat(varargin{3},size(cell2mat(cell_e(:)),1),1);
else
    connect = 0;
end

%preallocate
all_mtx = cell2mat(cell_e(:));
if ~exist('colors', 'var')
    try
        colors = linspecer(size(all_mtx,1));
    catch
        colors = repmat(.8.*[1 1 1], size(all_mtx,1), 1);
    end
end
all_mtx_grps = [];

means = nan(length(cell_e),1);
stds = nan(length(cell_e),1);
sqrt_l = nan(length(cell_e),1);

hold on;



cell_e_xaxis = cell(size(cell_e));
color_ct = 0;
for imtx = 1: length(cell_e)
    cell_e_xaxis{imtx} = nan(size(cell_e{imtx}));
    
    

        %colors = distinguishable_colors(length(cell_e{imtx}));
        mean_imtx = nanmean(cell_e{imtx});
        std_imtx = max([nanstd(cell_e{imtx}) realmin]);
        for idpt = 1:length(cell_e{imtx})
           
            %base_jitter_multiplier = length(cell_e{imtx})*0.0006 + 0.10;
            
            base_jitter_multiplier = length(cell_e{imtx})*0.0006 + 0.350;
            base_jitter = (rand(1)-0.5)*base_jitter_multiplier;
            dist_from_mean = abs(cell_e{imtx}(idpt) - mean_imtx);
            std_from_mean = dist_from_mean/std_imtx;
            bulb_correction = std_from_mean/7.5 + (rand(1)-0.5)*0.4;
            cell_e_xaxis{imtx}(idpt) = imtx+base_jitter*(1-bulb_correction);
            
            if exist('xaxis', 'var') && ~isempty(xaxis)
                %before = cell_e_xaxis{imtx}(idpt)
                cell_e_xaxis{imtx}(idpt) = cell_e_xaxis{imtx}(idpt) + diff([cell_e_xaxis{imtx}(idpt); xaxis(imtx)]) + base_jitter*(1-bulb_correction);
                %after = cell_e_xaxis{imtx}(idpt)
            end
            
           
            color_ct = color_ct+1;
            plot(cell_e_xaxis{imtx}(idpt), cell_e{imtx}(idpt), 'o', 'color', colors(color_ct,:), 'markersize', 5)
        end
    %end
    all_mtx_grps = [all_mtx_grps; repmat(imtx, size(cell_e{imtx}(:)))];
    means(imtx) = nanmean(cell_e{imtx});
    stds(imtx) = nanstd(cell_e{imtx});
    sqrt_l(imtx) = sqrt(sum(~isnan(cell_e{imtx}))); 
end
if connect == 1
    
    try
    cell_e_xaxis_mtx = cell2mat(cell_e_xaxis);
    catch
    end
    
    cell_e_mtx = cell2mat(cell_e);
    
    color_ct = 0;
    for subj = 1:size(cell_e_mtx,1)

        color_ct = color_ct+1;
        
        plot(cell_e_xaxis_mtx(subj,:), cell_e_mtx(subj,:), '-o', 'markersize', 5, 'color', colors(color_ct,:))  

    end
end
std_es = stds./sqrt_l;


set(gca,'TickLength',[0, 0]); box off
xlim([.5 length(cell_e)+.5])
%xticks(1:4)


if length(cell_e)==2
   if connect==1
       [~, pval, ~, stat_struct] = ttest(cell_e{1}, cell_e{2});
       title(['t(' num2str(stat_struct.df) ')=' num2str(round(abs(stat_struct.tstat).*1000)/1000) ', p=' num2str(round(pval.*1000)/1000)])
   else
       [~, pval, ~, stat_struct] = ttest2(cell_e{1}, cell_e{2});
       title(['t(' num2str(stat_struct.df) ')=' num2str(round(abs(stat_struct.tstat).*1000)/1000) ', p=' num2str(round(pval.*1000)/1000)])
   end
end


end

