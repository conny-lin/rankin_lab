function [SumData,legend] = get_speeddata_fromTrinity_v2(Trn,varargin)
%%
outputstyle = 'table';
vararginProcessor;
% process trinity file 
legend = {'wormid','time','tap','speed_bias'};

% find out total rows
[nrow, ~] = cellfun(@size,Trn(:,2));
nrowsum = sum(nrow);
SumData = nan(nrowsum,numel(legend));
rowstart = [0;cumsum(nrow)];
% get only time(1), tap(7), speed (4), bias (6)
wormid = Trn(:,1);
wormid = cellfun(@str2num,wormid);
for wrmi = 1:numel(wormid)
    wormid_current = wormid(wrmi);
    Data = Trn{wrmi,2};
    speed = Data(:,4);
    bias = Data(:,6);
    speed_dir = speed.*bias;
    time = Data(:,1);
    tap = Data(:,7);
    wid = repmat(wormid_current,size(Data,1),1);
%     pid = repmat(mwti,size(Data,1),1);        
%     D = [pid wid time tap speed_dir];
    D = [wid time tap speed_dir];

    SumData(rowstart(wrmi)+1:rowstart(wrmi+1),1:numel(legend)) = D;
end

switch outputstyle
    case 'table'
        SumData = array2table(SumData,'VariableNames',legend);
end
% SumData = array2table(SumData,'VariableNames',{'mwtid','wormid','time','tap','speed_bias'});