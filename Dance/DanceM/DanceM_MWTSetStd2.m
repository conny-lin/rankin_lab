function MWTSet = DanceM_MWTSetStd2(pMWT,funName,GroupVar,VarIn,pSave,outputsuffix)
% DanceM_MWTSetStd2 - use DanceM_createIndex2 for flexibility
% create MWTSet
MWTSet.ProgramName = funName;
MWTSet.timestamp = generatetimestamp; % generate timestamp
MWTSet.PATHS.pMWT = pMWT;

% record limits
if numel(VarIn) > 1
    callName = VarIn(1:2:numel(VarIn));
    setValue = VarIn(2:2:numel(VarIn));
    for x = 1:numel(callName)
        MWTSet.limit.(callName{x}) = setValue{x}; 
    end
end

% save folder
if isempty(outputsuffix) == 0
    pSave = [pSave,'/',funName,' ',outputsuffix];
else
    pSave = [pSave,'/',funName];
end
if isdir(pSave) == 0; mkdir(pSave); end
MWTSet.PATHS.pSaveA = pSave;

% create index
[MWTDB,MWTDBInd,VInd] = DanceM_createIndex2(pMWT); % create index


%% PROCESS VAR
% if VarName > 1
if numel(GroupVar) > 1
    [VInd,MWTDBInd,~,GroupVarName] = make_varref(MWTDBInd,GroupVar,VInd);
    GroupVarName = char(GroupVarName);
else
    GroupVarName = char(GroupVar);
end

% store variable in MWTSet
MWTSet.Info.MWTDB = MWTDB;
MWTSet.Info.MWTDBInd = MWTDBInd;
MWTSet.Info.VInd = VInd;
MWTSet.Info.GroupVarName = GroupVarName;


% clear temp file 
clearvars funName funPath VarIn;