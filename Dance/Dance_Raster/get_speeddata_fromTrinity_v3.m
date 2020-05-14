function [SumData,legend] = get_speeddata_fromTrinity_v3(Trn,varargin)
%%
outputstyle = 'table';
biascreen = 'exclude nan';
vararginProcessor;
% process trinity file 
legend = {'wormid','time','tap','bias','speed','speed_bias'};

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
    speed_dir = speed.*bias; % if bias is nan, no data
    time = Data(:,1);
    tap = Data(:,7);
    wid = repmat(wormid_current,size(Data,1),1);
    D = [wid time tap bias speed speed_dir];
    SumData(rowstart(wrmi)+1:rowstart(wrmi+1),1:numel(legend)) = D;
    
end

SumData = array2table(SumData,'VariableNames',legend);
SumData(isnan(SumData.bias),:) = [];
SumData(isnan(SumData.speed),:) = [];



