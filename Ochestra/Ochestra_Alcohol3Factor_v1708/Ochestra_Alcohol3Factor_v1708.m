function OcA = Ochestra_Alcohol3Factor_v1708(pMWT,pSave, varargin)
%% INITIALIZING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% clear; close all; 
% addpath('/Users/connylin/Dropbox/Code/Matlab/Library/General');
% pM = setup_std(mfilename('fullpath'),'RL','genSave',true); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% VARARGIN  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%20170815
% settings --------------------------------------------------------20170815
% control = 'N2';
% condition = {'NA','400mM'};
pDataBase = '/Volumes/COBOLT/MWT';
strain = '';
wildtype = 'N2';
overwrite = true;
pData = '';

%------------------------------------------------------------------20170815
% VARARGIN PROCESSOR  ---------------------------------------------20170815
vararginProcessor;
%------------------------------------------------------------------20170815
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%20170815


%% GET VARIABLES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%20170815
A = load(fullfile(pDataBase,'MWTDB.mat'));
if ~isempty(strain)
    B = A.MWTDB.strain;
    i = ismember(B.strain,strain);
    genotype = char(B.genotype(i)); % get genotype
    OcA.Info.genotype = genotype;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%20170815


%% DANCE ANALYSIS COMPONENTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%20170815
ci = true;
ri = true;
ai = true;
if overwrite==false && ~isempty(pData)
    p = fullfile(pData,'Dance_InitialEtohSensitivityPct_v1707','curve pct.csv');
    if exist(p,'file'); ci = false; end
    p = fullfile(pData,'Dance_ShaneSpark_v1707','Descriptive HabCurve.csv');
    if exist(p,'file'); ri = false; end
    p = fullfile(pData,'Dance_rType2_v1707','AccProb RMANOVA.txt');
    if exist(p,'file'); ai = false; end
end
if ci; Dance_InitialEtohSensitivityPct_v1707(pMWT,pSaveO); end
if ri; Dance_ShaneSpark_v1707(pMWT,pSaveO);end
if ai; Dance_rType2_v1707(pMWT,pSaveO);end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%20170815


%% GET GRAPH NUMBERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%20170815
Summary = struct;
% curve -----------------------------------------------------------20170815
varname = 'curve';
p = fullfile(pData,'Dance_InitialEtohSensitivityPct_v1707','curve pct.csv');
T = readtable(p);
% process group name
T.group = regexprep(T.group,'_400mM','');
% sort data
T = sortN2first(T.group,T);
% get data
A = struct;
A.G = T.group;
A.N = T.n;
A.Y = T.mean.*100;
A.E = T.se.*100;
% store in summary
Summary.(varname) = A;
%------------------------------------------------------------------20170815
% reversal --------------------------------------------------------20170815
varname = 'reversal';
p = fullfile(pData,'Dance_ShaneSpark_v1707','Dance_ShaneSpark_v1707.mat');
load(p);
% transform data
T = MWTSet.Graph.ByGroupPerPlate.RevFreq;
gori = T.groupname';
% get data
A = struct;
A.G = sortN2first(gori,T.groupname');
A.N = sortN2first(gori,T.N);
A.X = sortN2first(gori,T.tap);
A.Y = sortN2first(gori,T.Mean);
A.E = sortN2first(gori,T.SE);
% store in summary
Summary.(varname) = A;
%------------------------------------------------------------------20170815
% acceleration ----------------------------------------------------20170815
varname = 'acc';
p = fullfile(pData,'Dance_rType2_v1707','data.mat');
load(p);
% transform data
T = MWTSet.Graph.AccProb;
gori = T.groupname;
% get data
A = struct;
A.G = sortN2first(gori,T.groupname);
A.N = sortN2first(gori,T.N);
A.X = sortN2first(gori,T.X);
A.Y = sortN2first(gori,T.Y);
A.E = sortN2first(gori,T.E);
% store in summary
Summary.(varname) = A;
%------------------------------------------------------------------20170815
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%20170815


%% GRAPH %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%20170815
% setting  --------------------------------------------------------20170815
markersize = 3;
% gpos = {[0.07 0.11 0.1 0.45],...
%          [0.25 0.11 0.2 0.45],...
%          [0.55 0.11 0.2 0.45]};
% legpos = [0.04 0.66 0.1 0.22];
% textpos = [0.02 0.95 0.08 0.08];
w = 6;
h = 2.5;
gpos = {[0.07 0.05 0.14 0.85],...
         [0.31 0.05 0.28 0.85],...
         [0.71 0.05 0.28 0.85]};
     
legpos = [0.01 0.11 0.1 0.22];
textpos = [0.00 0.99 0.08 0.08];

close all
figure1 = figure('Visible','off');
%------------------------------------------------------------------20170815

% curve -----------------------------------------------------------20170815
% get data
varname = 'curve';
G = Summary.(varname);
% bar graph 
GN = G.G;
if ~strcmp(GN{1},'N2'); error('bad n2');end
Y = G.Y;
E = G.E;
% Create axes
axes1 = axes('Parent',figure1,'Position',gpos{1},'Box','off');
hold(axes1,'on');
% axes properties
set(axes1,'XTick',[1 2],'XTickLabel',{'+','-'});
xlim([0.5 2.5]);
ylim([-60 20]);
ylabel('% curve (etoh/control)')
bar1 = bar(Y,'Parent',axes1,'EdgeColor',[0 0.45 0.74],'FaceColor',[0 0.45 0.74]);
% errorbar
e1 = errorbar(Y,E,'LineStyle','none','Color','k','LineWidth',1,'Marker','none');
%------------------------------------------------------------------20170815


% reversal --------------------------------------------------------20170815
% get data
varname = 'reversal';
G = Summary.(varname);
GN = G.G;
X = G.X;
Y = G.Y;
E = G.E;
% Create axes
axes2 = axes('Parent',figure1,'Position',gpos{2});
hold(axes2,'on');
xlim([0 30.5]);
ylim([0,1]);
ylabel('P (reversal)');
xlabel('Tap');
e2 = errorbar(X,Y,E,'Marker','o','MarkerSize',markersize);
e2 = graphApplySetting(e2,'cathyline');
% Create multiple error bars using matrix input to errorbar
legname = {'+ 0mM','+ 400mM','(-) 0mM','(-) 400mM'};
for ei = 1:numel(GN)
    set(e2(ei),'DisplayName',legname{ei});
end
%------------------------------------------------------------------20170815

% Create legend %--------------------------------------------------20170815
% legend1 = legend(axes2,'show');
% set(legend1,'Position',legpos,'EdgeColor',[1 1 1]);
annotation(figure1,'textbox',textpos,...
    'String',{genotype},'FontWeight','bold','EdgeColor','none');
%------------------------------------------------------------------20170815


% acc -------------------------------------------------------------20170815
% get data
varname = 'acc';
G = Summary.(varname);
GN = G.G;
X = G.X;
Y = G.Y;
E = G.E;  
% Create axes
axes3 = axes('Parent',figure1,'Position',gpos{3});
hold(axes3,'on');
xlim([0 30.5])
ylim([0,1]);
ylabel('P (acceleration)');
xlabel('Tap');
e3 = errorbar(X,Y,E,'Marker','o','MarkerSize',markersize);
e3 = graphApplySetting(e3,'cathyline');
%------------------------------------------------------------------20170815

% save ------------------------------------------------------------20170815
savename = sprintf('%s %s',strain,genotype);

printfig(savename,pSave,'w',w,'h',h,'closefig',1);
%------------------------------------------------------------------20170815
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%20170815




%% STATS REPORT
% report_All(pData_Sen,pData_TAR,pData_SS,pSave,MWTDB,strain);
% ---------------

