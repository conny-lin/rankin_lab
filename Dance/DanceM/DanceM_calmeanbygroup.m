function [Result,gnamelist] = DanceM_calmeanbygroup(gname,Data)

gnamelist = unique(gname);
nRow = numel(gnamelist);
nCol = size(Data,2);
Result = nan(nRow,nCol);
for i = 1:nRow
    gn = gnamelist{i};
    gi = ismember(gname,gn);
    D = Data(gi,:);
    if size(D,1) > 1
       Result(i,:) = mean(D);
    else
        Result(i,:) = D;
    end
end