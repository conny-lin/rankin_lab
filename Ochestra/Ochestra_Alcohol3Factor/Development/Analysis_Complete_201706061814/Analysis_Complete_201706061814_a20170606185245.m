%% INITIALIZING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; close all; 
addpath('/Users/connylin/Dropbox/Code/Matlab/Library/General');
pM = setup_std(mfilename('fullpath'),'RL','genSave',true); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Function setting %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
strain = 'KP1097';
control = 'N2';
condition = {'NA','400mM'};
expdatemin = 20170210;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% DEFAULT SETTINGS  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% settings ++++++
pvsig = 0.05;
pvlim = 0.001;
w = 9;
h = 4;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Search MWTDB for experiments to analyze %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load(fullfile('/Volumes/COBOLT/MWT','MWTDB.mat')); % load MWTDB
MWTDBM = MWTDB.text; % pull MWTDB to new 
MWTDBM(MWTDBM.exp_date < expdatemin,:) = []; % remove older exp
i = ismember(MWTDBM.strain,strain); % get mwt plate with strain name
expnames = unique(MWTDBM.expname(i)); % get exp with mwt plate with strain name
MWTDBM(~ismember(MWTDBM.expname,expnames),:) = []; % remove exp without strain name
MWTDBM(~ismember(MWTDBM.strain,[control,strain]),:) = []; % retain only control strain
MWTDBM(~ismember(MWTDBM.rx,condition),:) = []; % retain only condition specified
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% get variable information %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pMWT = MWTDBM.mwtpath; % get mwtpath
genotype = MWTDB.strain.genotype(ismember(MWTDB.strain.strain,strain)); % get genotype
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% REVERSAL    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pSave_TWR = create_savefolder(fullfile(fileparts(pM),'Data',strain,'TWR')); % create save folder
MWTSet = Dance_ShaneSpark5(pMWT,pSave_TWR);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
return

%% #1: INITIAL SENSITIVITY %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pC = fullfile(pData_Strain,strain,'Etoh sensitivity','InitialEtohSensitivityPct','data.mat');
DC = load(pC,'DataMeta','MWTDB');
CS = CurveStats;
CS.mwtid = DC.DataMeta.mwtid;
CS.curve = DC.DataMeta.curve;
CS.MWTDB = DC.MWTDB;
[anovatxt,T] = anova(CS);
gn = regexprep(T.gnames,'_400mM','');


Summary.Curve.groupname = sortN2first(gn,gn);
Summary.Curve.N = sortN2first(gn,T.N);
Summary.Curve.Y = (sortN2first(gn,T.mean).*100)-100;
Summary.Curve.E = sortN2first(gn,T.SE).*100;
% --------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% ANALYSIS     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% chor ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% [Legend,pMWTpass,pMWTfailed] = chormaster5('ShaneSpark',pMWT); % chor
% fprintf('*** converting trinity.dat to .mat ***\n');
% convertTrinityDat2Mat(pMWTpass,0); % convert trinity to mat, and use only passed pMWT
% MWTDBLocal(ismember(MWTDBLocal.mwtpath,pMWTfailed),:) = []; % update MWTDB to exclude un-chor-able files
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

