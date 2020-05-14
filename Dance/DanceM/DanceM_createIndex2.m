function [MWTDb,MWTDbInd, VarIndex] = DanceM_createIndex2(pMWT)
%% MWTSet = createIndex(MWTSet) CREATE INDEX

%% MWTSet.Info.MWTDb 
% pMWT = MWTSet.PATHS.pMWT;
MWTDb = parseMWTinfo(pMWT);

%% MWTSet.Info.VarIndex 
% create index
a = fieldnames(MWTDb);
a(ismember(a,'Properties')) = [];
VarIndex = struct;
for x = 1:numel(a)
    VarIndex.(a{x}) = unique(MWTDb.(a{x}));
end
% MWTSet.Info.VarIndex = VarIndex;

%% MWTSet.Info.MWTDbInd 
% create index database
a = fieldnames(MWTDb);
a(ismember(a,'Properties')) = [];
fnames = a;
d = array2table(nan(size(MWTDb)),'VariableNames',fnames);
for x =1:numel(fnames)
    [~,j] = ismember(MWTDb.(fnames{x}), VarIndex.(fnames{x}));
    d.(fnames{x})= j;
end
MWTDbInd = d;
% MWTSet.Info.MWTDbInd = d;
