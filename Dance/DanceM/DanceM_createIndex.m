function MWTSet = DanceM_createIndex(MWTSet)
%% MWTSet = createIndex(MWTSet) CREATE INDEX

%% MWTSet.Info.MWTDb 
pMWT = MWTSet.PATHS.pMWT;
[p,fmwt] = cellfun(@fileparts,pMWT,'UniformOutput',0);
[p,gn] = cellfun(@fileparts,p,'UniformOutput',0);
[~,en] = cellfun(@fileparts,p,'UniformOutput',0);
MWTDb = table;
MWTDb.expname = en;
MWTDb.groupname = gn;
a = regexpcellout(gn,'_','split');
a = a(:,1);
MWTDb.strain = a;
a = regexpcellout(gn,'(?<=[_])\w*','match');
MWTDb.condition = a;
MWTDb.mwtname = fmwt;
MWTSet.Info.MWTDb = MWTDb;

%% MWTSet.Info.VarIndex 
% create index
a = fieldnames(MWTDb);
a(ismember(a,'Properties')) = [];
VarIndex = struct;
for x = 1:numel(a)
    VarIndex.(a{x}) = unique(MWTDb.(a{x}));
end
MWTSet.Info.VarIndex = VarIndex;

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
MWTSet.Info.MWTDbInd = d;
