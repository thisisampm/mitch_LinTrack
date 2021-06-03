function [cor_mtx, ebp_reference, ebp_adjacent] = sesh_to_sesh_corrs(tuning_curve_matrix, reference_session)
% compute average cell to cell correlation across sessions
%
% plots a correlation matrix and an errorbar plot 
%
% tuning_curve_matrix is a cell by bin by session rate matrix output by
% all_cell_tuning or common_cell_tuning functions. rows are same cells
% across sessions. a row can be nans if cell is inactive in that session.
%


% preallocate
cor_mtx = nan(size(tuning_curve_matrix,3)-1);
ebp_reference = nan(size(tuning_curve_matrix,1), size(tuning_curve_matrix,3)-1);
ebp_adjacent = nan(size(tuning_curve_matrix,1), size(tuning_curve_matrix,3)-1);

% for each pair of sessions
ref_sesh_ct = 0;
for isesh1 = 1:size(tuning_curve_matrix,3)
    for isesh2 = 1:size(tuning_curve_matrix,3)
                
        % for each cell
        cell_corrs = nan(size(tuning_curve_matrix,1),1);
        for ic = 1:size(tuning_curve_matrix,1)
            cell_corrs(ic) = corr(tuning_curve_matrix(ic,:,isesh1)', tuning_curve_matrix(ic,:,isesh2)');
        end
        
        % z transform rvals
        cell_corrs = atanh(cell_corrs);
        
        % reference session error bar plot info
        %
            % if adjacent session
            if isesh1==reference_session && isesh2>isesh1
                ref_sesh_ct = ref_sesh_ct +1;
                ebp_reference(:,ref_sesh_ct) = cell_corrs;
            end
            
        % adjacent error bar plot info
        %
            % if adjacent session
            if isesh2==(isesh1+1)
                ebp_adjacent(:,isesh1) = cell_corrs;
            end

        % correlation matrix info
        cor_mtx(isesh1, isesh2) = nanmean(cell_corrs);
        
    end
end

% plot corm
figure; 
imagesc(cor_mtx)
set(gca,'TickLength',[0, 0]); box off;
title('Average single cell correlations')
xlabel('Session')
ylabel('Session')


% plot errorbar of correlations to reference session
ebp_in = cell(1,size(ebp_reference,2));
for icomp = 1:size(ebp_reference,2)
    ebp_in{icomp} = ebp_reference(:,icomp);
end
figure; 
errorbar_plot(ebp_in, 0, [], .3.*[1 1 1]);
set(gca,'TickLength',[0, 0]); box off;
xticklabels(setdiff(1:size(tuning_curve_matrix,3), reference_session))
title(['Correlation with session ' num2str(reference_session)])
xlabel('Session')

[a b c d] = ttest2(ebp_in{1}, ebp_in{end})


% plot errorbar of adjacent sessions
ebp_in = cell(1,size(ebp_adjacent,2));
for icomp = 1:size(ebp_adjacent,2)
    ebp_in{icomp} = ebp_adjacent(:,icomp);
end
figure; 
errorbar_plot(ebp_in, 0, [], .3.*[1 1 1]);
set(gca,'TickLength',[0, 0]); box off;
title('Adjacent sessions')
xlabel('Comparison')


