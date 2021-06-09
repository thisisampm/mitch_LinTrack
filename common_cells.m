%cell_regist_mtx = cell_registered_struct.cell_to_index_map;
b = nan(size(cell_regist_mtx,2));
for i1 = 1:size(cell_regist_mtx,2)
    for i2 = 1:size(cell_regist_mtx,2)
        common_cell_sum = sum(cell_regist_mtx(:,i1)>0 & cell_regist_mtx(:,i2)>0);
        uncommon_cell_sum = sum([sum(cell_regist_mtx(:,i1)>0 & cell_regist_mtx(:,i2)==0) sum(cell_regist_mtx(:,i1)==0 & cell_regist_mtx(:,i2)>0)]);
        b(i1,i2) = common_cell_sum/(common_cell_sum + uncommon_cell_sum);
    end
end

figure; ampm_pcolor(b); title('Percentage of cells common to both sessions'); axis square; colorbar