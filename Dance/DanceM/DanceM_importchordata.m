function [MWTSet,Data] = DanceM_importchordata(MWTSet,varargin)

%% default
timerange = [];
chorname = fieldnames(MWTSet.legend_chor);
dancename = MWTSet.ProgramName;

%% vararginprocessor
vararginProcessor

%% post 
if isempty(timerange) && numel(timerange)==2 && timerange(1)<timerange(2)
    timerange_procses = false;
else
    timerange_procses = true;
end
%% input dealing 
if nargin<2
end

if ~iscell(chorname)
    chorname = {chorname};
end



%% process
pMWT = MWTSet.pMWT;

switch dancename
    case 'Dance_RapidTolerance_v1707'
        [Data,~,Legend2] = import_drunkposture2_dat(pMWT,'array');
        if numel(pMWT) ~=numel(Data); error('import number incorrect'); end
        Data = array2table(cell2mat(Data),'VariableNames',Legend2);
        % remove data outside of assay time
        if timerange_procses
            Data(Data.time < timerange(1) | Data.time > timerange(2),:) = [];
        end
        MWTSet.Raw = Data;
        MWTSet.legend_drunkposture2 = Legend2;
        
        
    otherwise
        error('no data import code available for %s',dancename);
end
