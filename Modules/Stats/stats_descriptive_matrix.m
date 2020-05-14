function A = stats_descriptive_matrix(B,col_array,col_label)


 % stats
n = sum(~isnan(B))';
mn = nanmean(B)';
mx = nanmax(B)';
mi = nanmin(B)';
sd = nanstd(B)';
% put in table
A = array2table([col_array n mn mi mx sd],...
    'VariableNames',{col_label,'n','mean','min','max','sd'});