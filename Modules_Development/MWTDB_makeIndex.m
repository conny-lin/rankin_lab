function MWTDB_makeIndex(Db)
%% MWTDB_makeIndex

%% FUNCTION PATHS
% general home path
pMatFun = '/Users/connylin/Dropbox/Code/Matlab';
pRLFun = '/Users/connylin/Dropbox/RL/Code/Modules';
% add packges
addpath([pMatFun,'/General']);
addpath([pRLFun,'/MWTDatabase']);
%% DEFAULTS 
nInput = 1;
%% VARARGIN PROCESSOR
if nargin > nInput
    A = varargin;
    if isinteger(numel(A)/2) == 1
        error('variable entries incorrect')
    end
    callName = A(1:2:numel(A));
    % check if call names are strings
    for x = 1:numel(callName)
       if ischar(callName{x}) == 0
           error('varargin incorrect');
       end
    end
    % get values
    setValue = A(2:2:numel(A));
    for x = 1:numel(callName)
        % determine input type
        if eval(sprintf('ischar(%s)',callName{x})) == 1
            eval(sprintf('%s = ''%s'';',callName{x},setValue{x}));
        elseif eval(sprintf('isnumeric(%s)',callName{x}))
            eval(sprintf('%s = %d;',callName{x},cell2mat(setValue(x))));
        elseif eval(sprintf('iscell(%s)',callName{x}))
            a = setValue{x};
            eval(sprintf('%s = a;',callName{x}));
        else
            error('need coding');
        end
    end
end


%% Code below not tested
return
%% GET PARSED 
Db.Properties.VariableNames

%% CREATE INDEX
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

% create index
a = fieldnames(MWTDb);
a(ismember(a,'Properties')) = [];
VarIndex = struct;
for x = 1:numel(a)
    VarIndex.(a{x}) = unique(MWTDb.(a{x}));
end
MWTSet.Info.VarIndex = VarIndex;

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