function [pMWT,MWTDB,Set] = search_MWTDB(varargin)


%% default
DBfname = 'MWTDB.mat';
pDB = '/Volumes/COBOLT/MWT';
ctrlexp = 'within';
ctrlgroup = {'N2'};

%% process
vararginProcessor_skipcheck;



%% GET MWT INFORMATION
% load database
load(fullfile(pDB, DBfname),'MWTDB');
% take out and rename database 
MWTDBM = MWTDB.text; 
% clear database file
clear MWTDB;
% get variable names from database
varnames = char(MWTDBM.Properties.VariableNames');
% process variables
Set = whos;
% take variable name out from input names
SetName = {Set(:).name};
% delete inputs do not match database variable names
Set(~ismember(SetName,varnames)) = [];
% get number of variables
setN = numel(Set);


%% QUARY DATABASE
VAL = true(size(MWTDBM,1),setN);
for seti = 1:setN
    v = char(Set(seti).name); % get variable name
    if ~eval(sprintf('isempty(%s);',v))
        str = sprintf('ismember(MWTDBM.%s, %s);',v,v);
        i = eval(str);
        VAL(~i,seti) = false;
    end
end
VAL2 = sum(VAL,2)==setN;
% sum(VAL2)
iMWTDB = VAL2;


%% ONLY INCLUDE CONTROL WITHIN EXPERIEMNTS
switch ctrlexp
    case 'within'
        % get group name
        g = unique(MWTDBM.groupname(iMWTDB));
        % get rid of control group name
        g(ismember(g,ctrlgroup)) = [];
        % index to group names not controls
        i = ismember(MWTDBM.groupname,g);
        % get expname
        e = unique(MWTDBM.expname(i));
        % get expname not matching 
        j = ~ismember(MWTDBM.expname,e);
        % exclude from quary
        iMWTDB(j) = false;
    otherwise
        
end


%% PREPARE OUTPUT
MWTDB = MWTDBM(iMWTDB,:);
pMWT = MWTDB.mwtpath;

end




















