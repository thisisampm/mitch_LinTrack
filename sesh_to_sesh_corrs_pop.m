function [cor_mtx, ebp_reference, ebp_adjacent] = sesh_to_sesh_corrs_pop(tuning_curve_matrix, reference_session)
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
ebp_reference = nan(size(tuning_curve_matrix,2), 1);
ebp_adjacent = nan(size(tuning_curve_matrix,2), 1);

% zscore
tuning_curve_matrix = zscore_mtx(tuning_curve_matrix);

% for each pair of sessions
ref_sesh_ct = 0;
for isesh1 = 1:size(tuning_curve_matrix,3)
    for isesh2 = 1:size(tuning_curve_matrix,3)
                
        % for each spatial bin
        bin_corrs = nan(size(tuning_curve_matrix,2),1);
        for ibin = 1:size(tuning_curve_matrix,2)
            rv1 = tuning_curve_matrix(:,ibin,isesh1);
            rv2 = tuning_curve_matrix(:,ibin,isesh2);
            nnan_idx = ~isnan(rv1) & ~isnan(rv2);
            bin_corrs(ibin) = corr(rv1(nnan_idx), rv2(nnan_idx));
        end
        
        % z transform rvals
        %bin_corrs_z = atanh(bin_corrs);
        
        % reference session error bar plot info
        %
            % if reference session
            if isesh1==reference_session && isesh2>isesh1
                ref_sesh_ct = ref_sesh_ct +1;
                ebp_reference(:,ref_sesh_ct) = bin_corrs;
            end
            
        % adjacent error bar plot info
        %
            % if adjacent session
            if isesh2==(isesh1+1)
                ebp_adjacent(:,isesh1) = bin_corrs;
            end

        % correlation matrix info
        cor_mtx(isesh1, isesh2) = nanmean(bin_corrs);
        
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


