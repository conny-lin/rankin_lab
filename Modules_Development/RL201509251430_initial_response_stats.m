%% RL201509251430_initial_response_stats
% find out which initial responses are significantly lower

%% get function path
addpath('/Users/connylin/Dropbox/MATLAB/Function_Library_Public');

pData = '/Users/connylin/Dropbox/Lab/Dance Output/20150925133008_Dance_ShaneSpark3';

%% get strains with significantly different initial
pFile = [pData,'/Stats Initial/RevFreq Initial posthoc bonferroni.csv'];
T = readtable(pFile);
prob = T.a3;
% find rows with no significant values
i = cellfun(@isempty,regexp(prob,'< 0.05'));
% keep only significant value table
T(i,:) = [];

% find two strains the same with each other
strain1 = T.a1;
a = regexpcellout(strain1,'_','split');
strain1_strain = a(:,1);
strain2 = T.a2;
a = regexpcellout(strain2,'_','split');
strain2_strain = a(:,1);
% index to rows for comparison of the same strain
i = cellfun(@strcmp,strain1_strain,strain2_strain);
% get strain with significantly different initial
strain_initial_diff = T.a1(i);


%% separate strains with lower or higher initial
% pFile = '/Users/connylin/Dropbox/Lab/Dance Output/20150925133008_Dance_ShaneSpark3/Tables/initial hab level summary.txt';
pFile = [pData,'/Graph Habituation curves/RevFreq.csv'];
T = readtable(pFile);
a = regexpcellout(T.groupname,' ','split');
strainlist = unique(a(:,1));
% put N2 first
strainlist(strcmp(strainlist,'N2')) = [];
strainlist = [{'N2'};strainlist];
a(cellfun(@isempty,a(:,2)),2) = {'0mM'};

D  = [];
for x = 1:numel(strainlist)
    strainname = strainlist{x};
    strainname_alcohol = [strainname ' 400mM'];
    a = T.t1(ismember(T.groupname, strainname));
    b = T.t1(ismember(T.groupname, strainname_alcohol));
    D(x,1) = b - a;  
end
A = table;
A.strainname = strainlist;
A.initial_diff = D;
A.sig = cell(size(A.strainname));
A.sig(ismember(A.strainname,strain_initial_diff)) = {'p < 0.05'};
A.sig(~ismember(A.strainname,strain_initial_diff)) = {'n.s.'};

%% save results
pSave = [pData,'/Stats Initial']; if isdir(pSave) == 0; mkdir(pSave); end
cd(pSave); 
writetable(A,'RevFreq Initial Report.csv');

%% separate graphs

