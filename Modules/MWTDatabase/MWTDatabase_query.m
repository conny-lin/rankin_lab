function DbT = MWTDatabase_query(pData,varargin)
%% MWTDatabase_query
% query MWTDatabase
% Input: 
%         pData = '/Volumes/COBOLT/MWT';
%     output settings
%         DisplayVar - variables asked to display on screen in cell arrays
%     Query terms
%         mwt_id   
%         mwt      
%         mwtpath  
%         expname  
%         exp_date 
%         tracker  
%         expter   
%         rc       
%         groupname = {'N2'};
%         strain   
%         rx       
%         preplate = 100;
%         ISI = 10;
%         postrec  
%         tapN = 30;
%         genotype
% 
%     Query selection terms
%         controlgroup = {'N2'};
%         gnameSearchType:
%             'any' (default) - groups can appear in any experiment
%             'within' - all groups mentioned must be within the same experiment 
%             'withcontrol' - groups must show up with controlgroup specified in controlgroup



%% FUNCTION PATHS
% general home path
% pMatFun = '/Users/connylin/Dropbox/Code/Matlab';
% pRLFun = '/Users/connylin/Dropbox/RL/Code/Modules';
% % add packges
% addpath([pMatFun,'/General']);
% addpath([pRLFun,'/MWTDatabase']);


%% DEFAULTS AND VARARGIN
nInput = 1;
% ISI = 10;
% tapN = 30;
% controlgroup = {'N2'};
DisplayVar = {'expname'};
% groupname = {'N2'};
% preplate = 100;
gnameSearchType = 'any';
% rc = '100s30x10s10s';

% varargin processer
% if nargin > nInput
%     A = varargin;
%     if isinteger(numel(A)/2) == 1
%         error('variable entries incorrect')
%     end
%     callName = A(1:2:numel(A));
%     % check if call names are strings
%     for x = 1:numel(callName)
%        if ischar(callName{x}) == 0
%            error('varargin incorrect');
%        end
%     end
%     % get values
%     setValue = A(2:2:numel(A));
%     for x = 1:numel(callName)
%         % determine input type
%         if eval(sprintf('ischar(%s)',callName{x})) == 1
%             eval(sprintf('%s = ''%s'';',callName{x},setValue{x}));
%         elseif eval(sprintf('isnumeric(%s)',callName{x}))
%             eval(sprintf('%s = %d;',callName{x},cell2mat(setValue(x))));
%         elseif eval(sprintf('iscell(%s)',callName{x}))
%             a = setValue{x};
%             eval(sprintf('%s = a;',callName{x}));
%         else
%             error('need coding');
%         end
%     end
% end
vararginProcessor;

%% LOAD DATABASE
D = load([pData,'/MWTDatabase.mat']);
MWTDatabase = D.MWTDatabase;
% get targets
Db = MWTDatabase.mwt;

%% PROCESS RESTRICTION VARIABLES
a = Db.Properties.VariableNames;
b = cell(size(a));
for x = 1:numel(a)
    if exist(a{x})== 1
        eval(sprintf('b{x} = %s;',a{x}));
    end
end
i = cellfun(@isempty,b);
b = [a;b];
b(:,i) = [];
Qvar = b;
fprintf('Query %d variables:\n',size(Qvar,2));
disp(char(Qvar(1,:)'))
fprintf('\n');
%% query
nVar = size(Qvar,2);
val = false(size(Db,1),nVar);
for x = 1:nVar
    varname = Qvar{1,x};
    terms = Qvar{2,x};
    val(:,x) = ismember(Db.(varname),terms);
end
i = sum(val,2) == nVar;
if sum(i) == 0;
   fprintf('No records matches query\n');
   return
end


%% process selection query
switch gnameSearchType
    case 'any'
        DbT = Db(i,:);
    case 'within'
        DbA = Db(i,:);
        eU = unique(DbA.expname);
        eUval = false(size(eU));
        for ei = 1:numel(eU)
            gu = unique(DbA.groupname(ismember(DbA.expname,eU{ei})));
            if (sum(ismember(groupname,gu)) == numel(groupname)) == true
                eUval(ei) = true;
            end
        end
        i = ismember(DbA.expname,eU(eUval)) & ...
            ismember(DbA.groupname,groupname);
        DbT = DbA(i,:);
    case 'withcontrol'
        DbA = Db(i,:);
        eU = unique(DbA.expname);
        eUval = false(size(eU));
        groupname = [groupname;controlgroup];
        for ei = 1:numel(eU)
            gu = unique(DbA.groupname(ismember(DbA.expname,eU{ei})));
            if (sum(ismember(groupname,gu)) == numel(groupname)) == true
                eUval(ei) = true;
            end
        end
        i = ismember(DbA.expname,eU(eUval)) & ...
            ismember(DbA.groupname,groupname);
        DbT = DbA(i,:);
end



%% create display output
A = DisplayVar;
for x = 1:numel(A)
    varname = A{x};
    varU = unique(DbT.(varname));
    disp(char(varU))
    fprintf('%s: %d\n',varname,numel(varU));
end

















