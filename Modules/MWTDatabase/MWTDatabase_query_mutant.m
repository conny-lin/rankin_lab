function [pMWT,DBT] = MWTDatabase_query_mutant(mutant,varargin)

%% mutant = cell array of mutant strains


%% default setting and process varargin
liquid_transfer = 1;
runcond = '100s30x10s10s';
cond = {'400mM'};
wildtype = 'N2'; 
vararginProcessor

%% get pMWT
% make groupname list
strain = [{wildtype} mutant];
nrow = numel(strain)*(numel(cond)+1);
groupnameList = cell(nrow,1);
n = 1;
for x = 1:numel(strain)
   groupnameList(n) = strain(x);
   n = n+1;
   groupnameList{n} = strjoin([strain(x),{'_'},cond],'');
   n = n+1;
end

% load database
load('/Users/connylin/Dropbox/rl/MWTDB/MWTDB.mat')
Database = MWTDB.text;
% get experiments containing mutant strains
i = ismember(Database.strain, mutant) ...
    & ismember(Database.rc, runcond);
% after liquid transfer protocol
if liquid_transfer
    i = i & Database.exp_date > 20111213;
end
% get experiment containing mutants
j = ismember(Database.expname, unique(Database.expname(i)));
DBT = Database(j,:);
% retain groups of interest
DBT(~ismember(DBT.groupname, groupnameList),:) = [];
pMWT = DBT.mwtpath; % get mwt paths
% reporting
fprintf('\ngroupname:\n');
disp(char(unique(DBT.groupname)))
fprintf('\nfrom exp:\n');
disp(char(unique(DBT.expname)))
fprintf('\n%d MWT folders\n',numel(pMWT));
