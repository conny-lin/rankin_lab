function [Result,ResultT] = DanceM_RawData_getInitial(Data)

msrlist = fieldnames(Data);
nRow = size(Data.(msrlist{1}),2);
nCol = numel(msrlist);
Result = nan(nRow,nCol);
for msri = 1:numel(msrlist)
    msrname = msrlist{msri};
    Result(:,msri) = Data.(msrname)(1,:)';
end
ResultT = array2table(Result,'VariableNames',msrlist);