function [mtx, sort_idx] = sort_rows_by_peak(mtx)
% sorts rows by the peak value
% earliest peak will be first row, latest peak will be last row
% if there are rows without variability, they will be added to bottom

% find peak vals
peak_pos = nan(size(mtx,1),1);
for irow = 1:size(mtx,1)
    if any(~isnan(mtx(irow,:)))
        peak_pos(irow) = find(mtx(irow,:)==max(mtx(irow,:)),1,'first');
    else
        peak_pos(irow) = size(mtx,2)+1;
    end
    if all(mtx(irow,:)==mtx(irow,1))
        peak_pos(irow) = size(mtx,2)+1;
    end
end

% sort
[~,sort_idx] = sort(peak_pos);
mtx = mtx(sort_idx,:);

