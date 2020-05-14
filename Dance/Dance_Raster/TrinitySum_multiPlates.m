function [plateSumm,recordTime,savename,MWTsourceID,D] = TrinitySum_multiPlates(pMWT,finish, start,varargin)

%% default
frameInt = 0.2;
saveoption =0;
savemasterfile = 0;
TrinityData = {};
displayopt = false;
vararginProcessor;

%% obtain trinity data
%% declare matrix
expectedFrameN = (finish-start)/frameInt;
plateSumm =[];
recordTime = [];
MWTsourceID = [];

D = table;
D.mwtpath = pMWT;
D.data = cell(numel(pMWT),1);

for x = 1:numel(pMWT)
    if displayopt
        processIntervalReporter(numel(pMWT),10,'-prcessing MWT files',x);
    end
    pMWTp = pMWT{x};
%     [~,mwtn] = fileparts(pMWTp);
    
    if isempty(TrinityData)
        [S,~,assayTime,d] = trinitySummary2(pMWTp,start,finish,...
            'frameInt',frameInt,'saveoption',saveoption,...
            'cleanup',0);
    else
        D = TrinityData{x,2};
        [S,~,assayTime,d] = trinitySummary2(pMWTp,start,finish,...
            'frameInt',frameInt,'saveoption',saveoption,...
            'cleanup',0,'D',D);
    
    end
    
    if ~isempty(S) % if has output
        % if directed to save master file
        if savemasterfile == 1; D.data{x} = d;end
        
        if size(S,2) > expectedFrameN % if more than expected frameN
%             warning('trinitySummary gave bigger frame N, delete extra [%s]',mwtn);
            S = S(:,1:expectedFrameN);
        elseif size(S,2) < expectedFrameN
%             warning('trinitySummary gave smaller frame N, skip plate [%s]',mwtn);
        else
            plateSumm = [plateSumm;S];
            n = size(S,1);
            MWTsourceID = [MWTsourceID;repmat(x,n,1)];
            recordTime = [recordTime;repmat(assayTime,n,1)];
        end
    end
end

%% save
% create save name
starttxt = num2str(start);
finishtxt =  num2str(round(finish));
ntext = num2str(size(plateSumm ,1));
savename = ['rasterPlot_',starttxt,'_',finishtxt,'_N_',ntext];
