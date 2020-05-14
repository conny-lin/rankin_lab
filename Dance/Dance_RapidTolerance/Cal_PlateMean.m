pMWT = MWTDB.mwtpath;
msr = Legend2(2:end);
x = Data.(Legend2{1});
% set up output
A = nan(numel(pMWT),numel(Legend2));
A(:,1) = MWTDB.mwtid;
for msri = 1:numel(msr)
    [gn,m] = grpstats(Data.(msr{msri}),x,{'gname','mean'});
    gn = cellfun(@str2num,gn);
    [i,j] = ismember(A(:,1),gn);
    A(i,msri+1) = m(j(i));
end
% transform to table
Sum = array2table(A,'VariableNames',Legend2);
% add table variables
S = innerjoin(MWTDB,Sum,'Keys','mwtid');
MWTSet.Data_Plate = Sum;
% export raw data
cd(pSave);
writetable(S,'raw.csv');