function MWTDB = makeMWTDatabase3(pData,pMWTDB)
%% makeMWTDatabase(pData,pMWTDB)
% INPUT:
%     pData = path to MWT data
%         pData = '/Volumes/COBOLT/MWT';
%     pMWTDB = path to MWT database storage site
%         pMWTDB = '/Volumes/COBOLT/MWT';

% run these scripts %%
% pData = '/Volumes/COBOLT/MWT';
% pMWTDB = '/Volumes/COBOLT/MWT';
% MWTDB = makeMWTDatabase3(pData,pMWTDB);


%% path
if nargin <2
   pMWTDB = pData; 
end

%% function path
addpath('/Users/connylin/Dropbox/Code/Matlab/Library/General');


%% load current MWTDB
cd(pMWTDB);
MWTDB_current = load('MWTDB.mat','MWTDB');


%% archive old MWTDB
% make new name
archivename = sprintf('MWTDB_%s',generatetimestamp);
% save old file as new name
eval(sprintf('%s = MWTDB_current.MWTDB;',archivename))
% archive
cd(pMWTDB);
eval(sprintf('save(''MWTDB.mat'',''%s'',''-append'')',archivename))


%% survey database
fprintf('\nGetting info from database %s... ',pData);
pMWTD = getMWTfile(pData);
fprintf('done\n');


%% GENERATE DATABASE -----------------------------------------------------
D = struct;
% MWT index
T = parseMWTinfo(pMWTD);
D.text = T;


%% group name
fprintf('getting groupname... ');
% unique group names
fGU = unique(D.text.groupname);
T= table;
T.groupname_id = [1:numel(fGU)]';
T.groupname = fGU;
T.expN = nan(numel(fGU),1);
T.mwtN = nan(numel(fGU),1);
for x = 1:numel(fGU)
    i = ismember(D.text.groupname,fGU{x}); 
    a = unique(D.text.expname(i));
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
    disp(char(unique(D.text.expname(ismember(D.text.groupname,fGproblem)))))
else
    fprintf('\ngroup names ok\n');
end


%% strain
fprintf('getting strain... ');
strain = D.text.strain;
strainU = unique(strain);
% validate
i = cellfun(@isempty,D.text.strain);
if sum(i) > 1
    fprintf('\nproblem strain names in experiment:\n');
    disp(char(unique(D.text.expname(i))));
    % display problem group names
    fprintf('\nproblem group names:\n');
    disp(char(unique(D.text.groupname(i))))
    % eliminate not qualified
    i = ~regexpcellout(strainU,'\<[A-Z]{1,}\d{1,}\>');
    strainU(i) = [];
else
    fprintf('\nstrain names ok\n');
end
fprintf('total %d valid strains\n',numel(strainU));
T = table;
T.strain_id = [1:numel(strainU)]';
T.strain = strainU;
D.strain = T;


%% rx
b = regexpcellout(D.text.rx,'_','split');
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
    i = ~cellfun(@isempty,regexp(D.text.groupname,['_',rxU{x}],'once'));
    a = unique(D.text.expname(i));
    T.expN(x) = numel(a);
    b = D.text.groupname(i);
    T.mwtN(x) = numel(b);
end
D.rxparts = T;


%% rc
fGU = unique(D.text.rc);
T = table;
T.rc_id = [1:numel(fGU)]';
T.rc = fGU;
T.expN = nan(numel(fGU),1);
T.mwtN = nan(numel(fGU),1);
for x = 1:numel(fGU)
    i = ismember(D.text.rc,fGU{x}); 
    a = unique(D.text.expname(i));
    T.expN(x) = numel(a);
    T.mwtN(x) = sum(i);
end
D.rc = T;


%% cond
rxU = unique(D.text.rx);
rxU = rxU(2:end);
T = table;
T.rx_id = [1:numel(rxU)]';
T.rx = rxU;
T.expN = nan(numel(rxU),1);
T.mwtN = nan(numel(rxU),1);
for x = 1:numel(rxU)
    i = ismember(D.text.rx,rxU{x}); 
    a = unique(D.text.expname(i));
    T.expN(x) = numel(a);
    T.mwtN(x) = sum(i);
end
D.rx = T;


%% exp report
a = cellfun(@strcat,D.text.rc, cellfunexpr(D.text.rc,'_'),D.text.groupname,'UniformOutput',0);
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
    i = ismember(D.text.rc,r{x}) & ismember(D.text.groupname,gn{x}); 
    a = unique(D.text.expname(i));
    T.expN(x) = numel(a);
    T.mwtN(x) = sum(i);
end
D.rc_groupname = T;


%% strain list - genotype
% SL = MWTDB_current.MWTDB.strainlist;
% get old strain/genotype list
SL = MWTDB_current.MWTDB.strain;
% 

%% FIX CODE FROM HERE
return
%%

D.strain = innerjoin(D.strain,SL);
A = D.strain;
[i,j] = ismember(D.text.strain,A.strain);
D.text.genotype = cell(size(D.text,1),1);
D.text.genotype(i) =  A.genotype(j(i));
if sum(~i) > 0
    fprintf('\n\nThese strains do not have genotype:\n');
    disp(unique(D.text.strain(~i)))
end


%% genotype
% find genes from genotype
a = D.strain.genotype;
b = regexp(a,'([a-z]+(-)\d+)|(wildtype*)|(?)','match');
% deal with exception (manual fix)
exception = {'C18H7.2'};
for x = 1:numel(exception)
    b(regexpcellout(a,exception{x})) = {exception(x)};
end
b(cellfun(@isempty,b)) = [];
% visulize what's not coded
ng = a(cellfun(@isempty,b)); 
disp(ng);
% process
b = celltakeout(b);
b = reshape(b,numel(b),1);
b(cellfun(@isempty,b)) = [];
b = unique(b);
b = [b;ng];
t = table;
t.gene_id = (1:numel(b))';
t.gene = b;
D.gene = t;


%% allele
a = D.strain.genotype;
b = regexp(a,'([a-z]+(-)\d+)|(wildtype*)|(?)','match');
% deal with exception
exception = {'C18H7.2'};
for x = 1:numel(exception)
    b(regexpcellout(a,exception{x})) = {exception(x)};
end
allele = cell(size(b));
for x = 1:size(b,1)
    c = b{x};
    allele{x} = cell(size(c));
    for y = 1:numel(c)
       if strcmp(c{y},'wildtype') == 1
           str = sprintf('(?<=(%s)\\()\\w+(?=\\))',c{y});
       else
          str = sprintf('(?<=(%s)\\()\\w+\\d+(?=\\))',c{y});
       end
       d = regexpcellout(a(x),str,'match');       
       allele{x}(y) = d;
    end
end
allele = celltakeout(allele);
% create allele list
ala = reshape(allele,numel(allele),1);
ala(cellfun(@isempty,ala)) = [];
t = table;
t.allele_id = (1:numel(ala))';
t.allele = ala;
D.allele = t;


%% allele - gene pair 
genes = b;
genes = celltakeout(genes);
a = cell(numel(genes),2);
for x = 1:numel(genes)
    a(x,1) = genes(x);
    a(x,2) = allele(x);
end
a(cellfun(@isempty,a(:,2)),:) = [];
% translate to id
genes = a(:,1);
allele = a(:,2);
i = ismember(genes,'?');
[i,j] = ismember(allele,D.allele.allele);
allele_id = j;
[i,j] = ismember(genes,D.gene.gene);
gene_id= j;
a = [gene_id allele_id];
a = unique(a,'rows');
t = table;
t.gene_id = gene_id;
t.allele_id = allele_id;
D.gene_allele_id = t;


%% gene sets
% exp set based on genes
% %% archive
a= {...
    'G protein',{'acy-1';'dgk-1';'egl-30';'egl-8';'goa-1';'gsa-1';'kin-2';'pde-4';'pkc-1';'ric-8'};
    'calcium',{'cnb-1'};
    'dopamine',{'cat-2', 'dop-1'};
    'synaptic',{'eat-16';'tom-1';'unc-13'};
    'nicotinic ion channel',{'eat-2'};
    'VG K channel',{'egl-2';'unc-103'};
    'neuropeptide Y',{'npr-1';'flp-18';'flp-20';'flp-21';'wildtyle;hawaiian'};
    'Glutamate',{'glr-1'};
    'neuroligin',{'mir-1';'nlg-1';'nrx-1';'pros-1'};
    'etoh insensitive',{'ptr-6'};
    'BK',{'slo-1'}
    };
D.expset = cell2table(a,'VariableNames',{'expset','genes'});


%% condition matrix


%% save
cd(pMWTDB);
MWTDB = D;
save('MWTDB.mat','MWTDB');
% writetable(D.rc_groupname,'rc_groupname.csv')
% writetable(D.text,'mwt.csv'); 



return


   
















