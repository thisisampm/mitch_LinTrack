function [cor_mtx_all, ebp_reference_all, ebp_adjacent_all] = sesh_to_sesh_corrs_multi_MDS(subj_cell, reference_session)
% compute average cell to cell correlation across sessions
%
% sesh_to_sesh_corrs but with values averaged within each subject before
% combining and plotting
%
% try loading subject_cell_ctl2wk.mat for inputs
%

% preallocate all subjs
%cor_mtx_all = nan(8, 8, length(subj_cell));
n_sesh = size(subj_cell{1},3); % Subj_cell enteries are a mouse array in cell x spatial bin x session matrix, so take the number of sessions from the last dimension of the first mouse 
cor_mtx_all = nan(n_sesh, n_sesh, length(subj_cell)); % Allocate a session x session x mouse matrix
ebp_reference_means = nan(length(subj_cell), n_sesh - 1);
ebp_adjacent_means = nan(length(subj_cell), n_sesh - 1);
ebp_reference_all = cell(length(subj_cell), n_sesh - 1);
ebp_adjacent_all = cell(length(subj_cell), n_sesh - 1);
    
% for each subject
for isubj = 1:length(subj_cell)

    
    % local subj matrix
    tuning_curve_matrix = subj_cell{isubj};
    
    % preallocate 1 subj
    cor_mtx = nan(size(tuning_curve_matrix,3)-1);
    ebp_reference = nan(size(tuning_curve_matrix,1), size(tuning_curve_matrix,3)-1);
    ebp_adjacent = nan(size(tuning_curve_matrix,1), size(tuning_curve_matrix,3)-1);

    % for each pair of sessions
    ref_sesh_ct = 0;
    for isesh1 = 1:size(tuning_curve_matrix,3)
        for isesh2 = 1:size(tuning_curve_matrix,3)

            % for each cell
            %{
            cell_corrs = nan(size(tuning_curve_matrix,1),1);
            for ic = 1:size(tuning_curve_matrix,1)
                cell_corrs(ic) = corr(tuning_curve_matrix(ic,:,isesh1)', tuning_curve_matrix(ic,:,isesh2)');
            end
            %}
            cell_corrs = diag(corr(tuning_curve_matrix(:,:,isesh1)',tuning_curve_matrix(:,:,isesh2)')); % corr will do a 1:1 Pearson correlation of all the columns in one matrix to the other, 
            % therefore the correlations between the same cells (same
            % position in each) will fall along the diagonal of the output.

            % reference session error bar plot info
            %
                % if ref session
                if isesh1==reference_session && isesh2>isesh1
                    ref_sesh_ct = ref_sesh_ct +1;
                    ebp_reference(:, ref_sesh_ct) = cell_corrs;
                    ebp_reference_all{isubj, ref_sesh_ct} = atanh(cell_corrs); % z transform
                end

            % adjacent error bar plot info
            %
                % if adjacent session
                if isesh2==(isesh1+1)
                    ebp_adjacent(:, isesh1) = cell_corrs;
                    ebp_adjacent_all{isubj, isesh1} = atanh(cell_corrs); % z transform
                end

            % correlation matrix info
            cor_mtx(isesh1, isesh2) = nanmean(cell_corrs);

        end
    end
    
    % load all subjs
    cor_mtx_all(:, :, isubj) = cor_mtx;
    ebp_reference_means(isubj,:) = nanmean(ebp_reference);
    ebp_adjacent_means(isubj,:) = nanmean(ebp_adjacent);
    
    

end

% plot corm
figure; 
imagesc(nanmean(cor_mtx_all,3))
set(gca,'TickLength',[0, 0]); box off;
title('Average subject single cell correlations')
xlabel('Session')
ylabel('Session')
axis square; colorbar


% SUBJECT AVERAGES
%
    % plot errorbar of correlations to reference session
    ebp_in = cell(1,size(ebp_reference_means,2));
    for icomp = 1:size(ebp_reference_means,2)
        ebp_in{icomp} = ebp_reference_means(:,icomp);
    end
    figure; 
    errorbar_plot([{ones(size(ebp_in{1}))} ebp_in], 1, []);
    set(gca,'TickLength',[0, 0]); box off;
    xticklabels(setdiff(1:8, reference_session))
    title(['Correlation with session ' num2str(reference_session) ' (subject means)'])
    xlabel('Session')

    % aesthetics
    hold on; plot(xlim, [0 0], 'k--')
    ylim([-1 1])

    [a b c d] = ttest2(ebp_in{1}, ebp_in{end})


    % plot errorbar of adjacent sessions
    ebp_in = cell(1,size(ebp_adjacent_means,2));
    for icomp = 1:size(ebp_adjacent_means,2)
        ebp_in{icomp} = ebp_adjacent_means(:,icomp);
    end
    figure; 
    errorbar_plot(ebp_in, 1, []);
    set(gca,'TickLength',[0, 0]); box off;
    title('Adjacent sessions (subject means)')
    xlabel('Comparison')

    % aesthetics
    hold on; plot(xlim, [0 0], 'k--')
    ylim([-1 1])
    
    
    
% ALL CELLS
%

    % plot errorbar of correlations to reference session
    figure; 
    errorbar_plot_multi(ebp_reference_all, 0, [], .3.*[1 1 1]);
    set(gca,'TickLength',[0, 0]); box off;
    xticklabels(setdiff(1:size(tuning_curve_matrix,3), reference_session))
    title(['Correlation with session ' num2str(reference_session) ' (all cells)'])
    xlabel('Session')

    % aesthetics
    hold on; plot(xlim, [0 0], 'k--')
    ylim([-1 1])
    maxr = .99; tic_vect = [-maxr -.9 -.8 -.6 -.3 0 .3 .6 .8 .9 maxr];
    ylim(atanh([-maxr maxr])); 
    ylim_hold = ylim;
    ylim([ylim_hold(1) ylim_hold(2)+(ylim_hold(2)*.25)])
    yticks(atanh(tic_vect)); yticklabels(tic_vect)

    [a b c d] = ttest2(ebp_in{1}, ebp_in{end})


    % plot errorbar of adjacent sessions
    figure; 
    errorbar_plot_multi(ebp_adjacent_all, 0, [], .3.*[1 1 1]);
    set(gca,'TickLength',[0, 0]); box off;
    title('Adjacent sessions (all cells)')
    xlabel('Comparison')

    % aesthetics
    hold on; plot(xlim, [0 0], 'k--')
    ylim([-1 1])
    maxr = .99; tic_vect = [-maxr -.9 -.8 -.6 -.3 0 .3 .6 .8 .9 maxr];
    ylim(atanh([-maxr maxr])); 
    ylim_hold = ylim;
    ylim([ylim_hold(1) ylim_hold(2)+(ylim_hold(2)*.25)])
    yticks(atanh(tic_vect)); yticklabels(tic_vect)




