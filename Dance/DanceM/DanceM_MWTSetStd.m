
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
MWTSet = DanceM_createIndex(MWTSet); % create index

% clear temp file 
clearvars funName funPath VarIn;