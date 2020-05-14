function D = get_speeddata_fromTrinity(Trinity)

% process trinity file 
SumDatap = cell(size(Trinity,1));
for mwti = 1:numel(Trinity,1);
    % find out total rows
    [nrow, ~] = cellfun(@size,Trinity{mwti}(:,2));
    nrowsum = sum(nrow);
    SumData = nan(nrowsum,5);
    rowstart = [0;cumsum(nrow)];
    % get only time(1), tap(7), speed (4), bias (6)
    wormid = Trinity{mwti}(:,1);
    wormid = cellfun(@str2num,wormid);
    for wrmi = 1:numel(wormid)
        wormid_current = wormid(wrmi);
        Data = Trinity{mwti}{wrmi,2};
        speed = Data(:,4);
        bias = Data(:,6);
        speed_dir = speed.*bias;
        time = Data(:,1);
        tap = Data(:,7);
        wid = repmat(wormid_current,size(Data,1),1);
        pid = repmat(mwti,size(Data,1),1);        
        D = [pid wid time tap speed_dir];
        SumData(rowstart(wrmi)+1:rowstart(wrmi+1),1:5) = D;
    end
    SumDatap{mwti} = SumData;
end
% put data in D
D = cell2mat(SumDatap);
clear Trinity SumDatap SumData;
% convert to table
D = array2table(D,'VariableNames',{'mwtid','wormid','time','tap','speed_bias'});