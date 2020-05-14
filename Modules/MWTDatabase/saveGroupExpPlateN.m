function saveGroupExpPlateN(MWTDB,pSave)


%%
gnu = unique(MWTDB.groupname);

T = table;
T.groupname = gnu;
T.expN = nan(size(T,1),1);
T.plateN = nan(size(T,1),1);
for gi = 1:numel(gnu)
   i = ismember(MWTDB.groupname,gnu(gi));
   T.expN(gi) = numel(unique(MWTDB.expname(i)));
   T.plateN(gi) = sum(i);
end
cd(pSave);
writetable(T,'group_N.csv');