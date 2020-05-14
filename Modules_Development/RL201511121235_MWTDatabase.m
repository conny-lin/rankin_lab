% load cobolt database generated from RL201511101653_tsfAnalysisResults

addpath('/Users/connylin/Dropbox/MATLAB/Function_Library_Public');
addpath('/Users/connylin/Dropbox/MATLAB/Programs_RankinLab/Library/Modules/MWTDatabase');

%% get data
pHome = '/Volumes/COBOLT/MWT';
[pMWTD] = getMWTfile(pHome);
[pG,fMWT] = cellfun(@fileparts,pMWTD,'UniformOutput',0);
[pE,fG] = cellfun(@fileparts,pG,'UniformOutput',0);
[~,fE] = cellfun(@fileparts,pE,'UniformOutput',0);

D = struct;
%% generate database
% MWT index
T = table;
T.mwt_id = [1:numel(fMWT)]';
T.mwt = fMWT;
T.mwtpath = regexprep(pMWTD,pHome,'');
T.expname = fE;
% T.expset = 
T.exp_date = cellfun(@str2num,regexpcellout(fE,'\<\d{8}','match'));
T.tracker = char(regexpcellout(fE,'(?<=\<\d{8})[A-Z]','match'));
T.expter = char(regexpcellout(fE,'(?<=\<\d{8}[A-Z][_])[A-Z]{2}','match'));
a = regexpcellout(fE,'_','split');
T.rc = a(:,3);
T.groupname = fG;
T.strain = regexpcellout(fG,'\<[A-Z]{1,}\d{1,}','match');
a = regexpcellout(fG,'_','split');
T.rx = regexpcellout(fG,'(?<=\<[A-Z]{1,}\d{1,}_)\w{1,}','match');
T.rx(cellfun(@isempty,T.rx)) = {'NA'};
a = regexpcellout(T.rc,'\d{1,}(?=s)','match');
T.preplate = cellfun(@str2num,a(:,1));
T.ISI = cellfun(@str2num,a(:,2));
T.postrec = cellfun(@str2num,a(:,3));
T.tapN = cellfun(@str2num,regexpcellout(T.rc,'\d{1,}(?=x)','match'));
D.mwt = T;



% group name
% unique group names
fGU = unique(D.mwt.groupname);
T= table;
T.groupname_id = [1:numel(fGU)]';
T.groupname = fGU;
T.expN = nan(numel(fGU),1);
T.mwtN = nan(numel(fGU),1);
for x = 1:numel(fGU)
    i = ismember(D.mwt.groupname,fGU{x}); 
    a = unique(D.mwt.expname(i));
    T.expN(x) = numel(a);
    T.mwtN(x) = sum(i);
end
T.problem = false(numel(fGU),1);
T.problem(~regexpcellout(fGU,'^[A-Z]{1,}\d{1,}') & ~regexpcellout(fGU,'\<Test\>')) = true;
D.groupname = T;
% report problem
if sum(T.problem == true) > 1
    fprintf('\nproblem group names:\n');
    fGproblem = T.groupname(T.problem == true);
    disp(char(fGproblem));
    fprintf('\nfrom exp:\n');
    disp(char(unique(D.mwt.expname(ismember(D.mwt.groupname,fGproblem)))))
else
    fprintf('\ngroup names ok\n');
end



% strain
strain = D.mwt.strain;
strainU = unique(strain);
% validate
i = cellfun(@isempty,D.mwt.strain);
if sum(i) > 1
    fprintf('\nproblem strain names:\n');
    disp(char(D.mwt.expname(i)));
    % eliminate not qualified
    i = regexpcellout(strainU,'\<[A-Z]{1,}\d{1,}\>');
    strainU(i) = [];
else
    fprintf('\nstrain names ok\n');
end
fprintf('total %d valid strains\n',numel(strainU));
T = table;
T.strain_id = [1:numel(strainU)]';
T.strain = strainU;
D.strain = T;



% rx
b = regexpcellout(D.mwt.rx,'_','split');
rx = {};
for x = 1:numel(b)
   rx = [rx;b{x}];
end
rxU = unique(rx);
T = table;
T.rxparts_id = [1:numel(rxU)]';
T.rxparts = rxU;
T.problem = false(numel(rxU),1);
rxProblemIgore = {'0A','0B','2A','2B','3A','3B'};
i = ismember(rxU,rxProblemIgore);
T.problem(i) = true;

T.expN = nan(numel(rxU),1);
T.mwtN = nan(numel(rxU),1);
T.solution = cell(numel(rxU),1);
for x = 1:numel(rxU)
    i = ~cellfun(@isempty,regexp(D.mwt.groupname,['_',rxU{x}],'once'));
    a = unique(D.mwt.expname(i));
    T.expN(x) = numel(a);
    b = D.mwt.groupname(i);
    T.mwtN(x) = numel(b);
end
D.rxparts = T;

% rc
fGU = unique(D.mwt.rc);
T = table;
T.rc_id = [1:numel(fGU)]';
T.rc = fGU;
T.expN = nan(numel(fGU),1);
T.mwtN = nan(numel(fGU),1);
for x = 1:numel(fGU)
    i = ismember(D.mwt.rc,fGU{x}); 
    a = unique(D.mwt.expname(i));
    T.expN(x) = numel(a);
    T.mwtN(x) = sum(i);
end
D.rc = T;

% cond
rxU = unique(D.mwt.rx);
rxU = rxU(2:end);
T = table;
T.rx_id = [1:numel(rxU)]';
T.rx = rxU;
T.expN = nan(numel(rxU),1);
T.mwtN = nan(numel(rxU),1);
for x = 1:numel(rxU)
    i = ismember(D.mwt.rx,rxU{x}); 
    a = unique(D.mwt.expname(i));
    T.expN(x) = numel(a);
    T.mwtN(x) = sum(i);
end
D.rx = T;

% exp report
a = cellfun(@strcat,D.mwt.rc, cellfunexpr(D.mwt.rc,'_'),D.mwt.groupname,'UniformOutput',0);
A = unique(a);
gn = regexpcellout(A,'(?<=[_])\w{1,}','match');
r = regexpcellout(A,'_','split');
r = r(:,1);



T = table;
T.rc_groupname_id = [1:numel(A)]';
T.rc_groupname = A;
T.rc = r;
T.groupname = gn;
T.expN = nan(numel(A),1);
T.mwtN = nan(numel(A),1);
for x = 1:numel(A)
    i = ismember(D.mwt.rc,r{x}) & ismember(D.mwt.groupname,gn{x}); 
    a = unique(D.mwt.expname(i));
    T.expN(x) = numel(a);
    T.mwtN(x) = sum(i);
end
D.rc_groupname = T;


% save
cd(pHome);
MWTDatabase = D;
save('MWTDatabase.mat','MWTDatabase');
writetable(D.rc_groupname,'rc_groupname.csv')
writetable(D.mwt,'mwt.csv');    

return






%% address exp naming problem

a = unique(D.mwt.groupname);
i = regexpcellout(a,'\<[A-Z]{1,}\d{1,}');
fGProblem = a(~i);

i = ismember(D.mwt.groupname,fGProblem);
unique(D.mwt.expname(i))
%%
for x = 5%:numel(fGProblem)
    fGProblem(x)
    i = ismember(D.mwt.groupname,fGProblem(x));
    unique(D.mwt.expname(i))
end


















