%% RL201509251430_initial_response_stats
% find out which initial responses are significantly lower

%% get function path
clear;
addpath('/Users/connylin/Dropbox/MATLAB/Function_Library_Public');
pData = '/Users/connylin/Dropbox/Lab/Dance Output/20150925133008_Dance_ShaneSpark3';

%% get strains with significantly different initial
pFile = [pData,'/matlab.mat'];
Data = load(pFile,'MWTSet');
D = Data.MWTSet.Data.Raw.Y.RevFreq;
%%
pMWT = Data.MWTSet.MWTInfo.pMWT;
[p,platename] = cellfun(@fileparts,pMWT,'UniformOutput',0);
[p,groupname] = cellfun(@fileparts,p,'UniformOutput',0);
[~,expname] = cellfun(@fileparts,p,'UniformOutput',0);
% combine group and exp name
groupname_expname = cellfun(@strcat,groupname,cellfunexpr(expname,'**'),expname,'UniformOutput',0);
% create strain name
a = regexpcellout(groupname,'_','split');
strainname = a(:,1);
a(:,2) = regexprep(a(:,2),'400MM','400mM');
% create alohol condition
i = cellfun(@isempty,a(:,2)); % index to alcohol = 0mM
a(i,2) = {'0mM'};
alcohol_condition = a(:,2);

%% organize by strain x alcohol condition
strainname_list = unique(strainname);
% put N2 first
strainname_list(ismember(strainname_list,'N2')) = [];
strainname_list = [{'N2'};strainname_list];
strainname_list_alcohol = cellfun(@strcat,...
        strainname_list,cellfunexpr(strainname_list,'_400mM'),...
        'UniformOutput',0);

%% organize data by condition
alcohol_condition_list = unique(alcohol_condition);
display 'alcohol condition:'; disp(char(alcohol_condition_list));

%% get initial
D_initial = D(1,:)';

%% caluclate mean for each exp/group by alcohol condition
expname_list = unique(expname);
rowN = 1;
A = struct;
for x = 1:numel(strainname_list)
    sn = strainname_list(x);
    sna = strainname_list_alcohol(x);
    fprintf('strain: %s\n',char(sn));
    % validate equalness
    a = regexpcellout(sna,'_','split');
    if strcmp(sn,a(:,1)) == 0; 
        error('strain name is not equal');
    end

    % get per exp
    
    for e = 1:numel(expname_list)
        en = expname_list(e);
        j = ismember(expname,en);

        % get 0mM
        s = sn;

    
        i = ismember(groupname,s);
        k = i == true & j == true;
        if sum(k) > 0
            A.strain(rowN,1) = s;
            A.expname(rowN,1) = en;
            d = D(1,k);
            A.N(rowN,1) = numel(d);
            A.mean(rowN,1) = mean(d);
            A.SE(rowN,1) = std(d)./sqrt((numel(d)-1));
        end
        
        i = ismember(groupname,sna);
        k = i == true & j == true;
        if sum(k) > 0
            d = D(1,i == true & j == true);
            A.N(rowN,2) = numel(d);
            A.mean(rowN,2) = mean(d);
            A.SE(rowN,2) = std(d)./sqrt((numel(d)-1));
        end
        if sum(k) > 0
            rowN = rowN + 1;
            fprintf('--expname: %s\n',char(en));
        end
    end
  
end

%% calculate mean for each exp/group
strainU = unique(A.strain);
%% assign color to each strain
n = numel(strainU);
% s = repmat(10,A.strain,1);
% load seamount;
% cl = z(1:n);
% c = zeros(n,1);
% for si = 1:numel(strainU)
%     strainname = strainU{si};
%     i = ismember(A.strain,strainname);
%     c(i) = cl(si);    
% end

%%
X = A.mean(:,1);
Y = A.mean(:,2);
g = A.strain;
fig = gscatter(X,Y,g);

pSave = [pData,'/Scatter TEMP']; if isdir(pSave) == 0; mkdir(pSave); end
cd(pSave);
savefigeps(mfilename,pSave)



return
%%
% c = linspace(1,10,length(X));
s = 30;
  
for si = 1%:numel(strainU)
    strainname = strainU{si};
    i = ismember(A.strain,strainname);
    X = A.mean(i,1);
    Y = A.mean(i,2);
    scatter(X,Y,s,[0 0 0],'fill')
end



% plot regression

% repeat for hab level
% repeat for hab rate 1
% repeat for hab rate 2




%%
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

