function Ochestra_Alcohol3Factor(strain,pDHome,varargin)
%% INITIALIZING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% clear; close all; 
% addpath('/Users/connylin/Dropbox/Code/Matlab/Library/General');
% pM = setup_std(mfilename('fullpath'),'RL','genSave',true); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Function setting %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% strain = 'KP1097';
% pDHome = '/Users/connylin/Dropbox/RL Data by strain';


%% VARARGIN  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
control = 'N2';
condition = {'NA','400mM'};
expdatemin = 20140210;
pSave = fullfile(pDHome,strain);
pDataBase = '/Volumes/COBOLT/MWT';
pvsig = 0.05;
pvlim = 0.001;
w = 9;
h = 4;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% VARARGIN PROCESSOR  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
vararginProcessor;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%% Search MWTDB for experiments to analyze %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load(fullfile(pDataBase,'MWTDB.mat')); % load MWTDB
MWTDBM = MWTDB;
MWTDB = MWTDBM.text; % pull MWTDB to new 
MWTDB(MWTDB.exp_date < expdatemin,:) = []; % remove older exp
i = ismember(MWTDB.strain,strain); % get mwt plate with strain name
expnames = unique(MWTDB.expname(i)); % get exp with mwt plate with strain name
MWTDB(~ismember(MWTDB.expname,expnames),:) = []; % remove exp without strain name
MWTDB(~ismember(MWTDB.strain,[{control},{strain}]),:) = []; % retain only control strain
MWTDB(~ismember(MWTDB.rx,condition),:) = []; % retain only condition specified
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% get variable information %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pMWT = MWTDB.mwtpath; % get mwtpath
genotype = char(MWTDBM.strain.genotype(ismember(MWTDBM.strain.strain,strain))); % get genotype
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% DANCE ANALYSIS COMPONENTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INITIAL SENSITIVITY 
Data_Sen = Dance_InitialEtohSensitivityPct(pMWT,'pSave',pSave); % run initial sensitivity
pData_Sen = Data_Sen.PATHS.pSaveA;
pC = fullfile(pData_Sen,'data.mat');
DC = load(pC,'DataMeta','MWTDB');
CS = CurveStats;
CS.mwtid = DC.DataMeta.mwtid;
CS.curve = DC.DataMeta.curve;
CS.MWTDB = DC.MWTDB;
[~,T] = anova(CS);
gn = regexprep(T.gnames,'_400mM','');

Summary.Curve.groupname = sortN2first(gn,gn);
Summary.Curve.N = sortN2first(gn,T.N);
Summary.Curve.Y = (sortN2first(gn,T.mean).*100)-100;
Summary.Curve.E = sortN2first(gn,T.SE).*100;
% --------------------------


% REVERSAL
Data_SS = Dance_ShaneSpark5(pMWT,'pSave',pSave);
pData_SS = Data_SS.PATHS.pSaveA;
load(fullfile(pData_SS,'Dance_ShaneSpark5.mat')); % load data
D = MWTSet.ByGroupPerPlate.RevFreq;
R = TWR_graph_data(D);
Summary.TWR = R;
% -------------------------------


% ACCELERATION
Data_TAR = Dance_rType2(pMWT,'pSave',pSave,'genotype',genotype); % rType
pData_TAR = Data_TAR.PATHS.pSaveA;
A = load(fullfile(pData_TAR,'data.mat'));
Summary.TAR = A.G.AccProb;
% -----------------------

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% GRAPH %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
markersize = 4.5;
gp = graphsetpack('cathyline');
gpn = fieldnames(gp);

close all
figure1 = figure('Visible','off');

% bar graph +++++++++++++
GN = Summary.Curve.groupname;
if ~strcmp(GN{1},'N2'); error('bad n2');end
Y = Summary.Curve.Y;
E = Summary.Curve.E;
% Create axes
axes1 = axes('Parent',figure1,'Position',[0.07 0.11 0.1 0.45],'Box','off');
hold(axes1,'on');
% axes properties
set(axes1,'XTick',[1 2],'XTickLabel',{'+','-'});
xlim([0.5 2.5]);
ylim([-60 20]);
ylabel('% curve (etoh/control)')
bar1 = bar(Y,'Parent',axes1,'EdgeColor',[0 0.447058826684952 0.74117648601532],'FaceColor',[0 0.447058826684952 0.74117648601532]);
% errorbar
e1 = errorbar(Y,E,'LineStyle','none','Color','k','LineWidth',1,'Marker','none');
% ------------------------


% plot 2 ++++++++++++++++
GN = Summary.TWR.gn;
X = Summary.TWR.X;
Y = Summary.TWR.Y;
E = Summary.TWR.E;
if ~strcmp(GN{1},'N2')
    X = sortN2first(GN,X')';
    Y = sortN2first(GN,Y')';
    E = sortN2first(GN,E')';
    GN = sortN2first(GN,GN);
end 
% Create axes
axes2 = axes('Parent',figure1,'Position',[0.25 0.11 0.3 0.8]);
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
% Create legend
legend1 = legend(axes2,'show');
set(legend1,'Position',[0.04 0.66 0.1 0.22],'EdgeColor',[1 1 1]);
annotation(figure1,'textbox',[0.02 0.95 0.08 0.08],...
    'String',{genotype},'FontWeight','bold','EdgeColor','none');
% ------------------------


% plot 3 +++++++++++
GN = Summary.TAR.groupname;
X = Summary.TAR.X;
Y = Summary.TAR.Y;
E = Summary.TAR.E;
if ~strcmp(GN{1},'N2')
    X = sortN2first(GN,X')';
    Y = sortN2first(GN,Y')';
    E = sortN2first(GN,E')';
    GN = sortN2first(GN,GN);
end    
% Create axes
axes3 = axes('Parent',figure1,'Position',[0.65 0.11 0.3 0.8]);
hold(axes3,'on');
xlim([0 30.5])
 ylim([0,1]);

ylabel('P (acceleration)');
xlabel('Tap');
e3 = errorbar(X,Y,E,'Marker','o','MarkerSize',markersize);
e3 = graphApplySetting(e3,'cathyline');
% ------------------------

% save -------------------------------------------
printfig(strain,pSave,'w',w,'h',h,'closefig',1);
% -------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%% STATS REPORT
% curve sensitivity +++++++++++++++
curveData = rmANOVACurveReport;
curveData.datapath = fullfile(pData_Sen,'curve MANOVA.txt');
% calculate percentage difference 
pC = fullfile(pData_Sen,'data.mat');
DC = load(pC,'DataMeta','MWTDB');
CS = CurveStats;
CS.mwtid = DC.DataMeta.mwtid;
CS.curve = DC.DataMeta.curve;
CS.MWTDB = DC.MWTDB;
txt = reportCurveSensitivity(curveData,CS);
% ---------------------------------

% TAR ++++++++++
TAR = rmANOVATARReport;
TAR.datapath = fullfile(pData_TAR, 'AccProb RMANOVA.txt');
txt2 = reportTAR(TAR);
% --------------

% TWR  ++++++++++
TWR = rmANOVATWR;
TWR.datapath = fullfile(pData_SS, 'RMANOVA.txt');
txt3 = reportTWR(TWR);
% --------------


txtFinal = sprintf('%s\n%s\n%s',txt,txt2,txt3);


% EXPORT REPORT +++++++++
p2 = fullfile(pSave,sprintf('%s stats writeup.txt',strain));
fid = fopen(p2,'w');
fprintf(fid,'%s',txtFinal);
fclose(fid);
%-----------

% save objects ++++
p3 = fullfile(pSave,sprintf('%s.mat',strain));
save(p3,'curveData','TAR','TWR','MWTDB');
% ---------------

