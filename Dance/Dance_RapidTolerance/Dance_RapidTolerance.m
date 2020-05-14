function MWTSet = Dance_RapidTolerance(pMWT,varargin)
%% Dance_DrunkMoves_RapidTolerance
% updated: 201603270930

%% Paths
addpath(fileparts(mfilename('fullpath')));

%% DEFAULTS ---------------------------------------------------------------
fprintf('\n** Running %s **\n',mfilename);
pSave = '/Users/connylin/Dropbox/RL/Dance Output';
%%
pDataBase = '/Volumes/COBOLT/MWT/MWTDB.mat';
timeStartSet = 240;
timeEndSet = 300;
% varargin
vararginProcessor;
% output
MWTSet.Input.pMWT = pMWT;
MWTSet.Input.timeStartSet = timeStartSet;
MWTSet.Input.timeEndSet = timeEndSet;
% create save
pSave = [pSave,'/',mfilename]; if isdir(pSave)==0; mkdir(pSave); end
% save a copy of this to pSave
copyfile([mfilename('fullpath'),'.m'], [pSave,'/',mfilename,'_a',generatetimestamp,'.m']);


%% PREPARE ----------------------------------------------------------------
% CHOR
[~,~,pfailed] = chormaster5('DrunkPostureOnly',pMWT);
% only pMWT with chor can be analyzed
if isempty(pfailed)==false
   warning('these mwt files do not have chor');
   a = parseMWTinfo(pfailed);
   disp(char(a.mwtname))
end
% update MWTDB from database
pMWT = pMWT(~ismember(pMWT,pfailed)); 
Db = load(pDataBase);
Db = Db.MWTDB.text;
i = ismember(Db.mwtpath,pMWT);
if sum(i) ~= numel(pMWT)
   error('database missing input pMWT info');
end
[i,j] = ismember(pMWT,Db.mwtpath);
MWTDB = Db(j(i),:); clear Db;
MWTDB.mwtid_db = MWTDB.mwtid;
MWTDB.mwtid = [1:size(MWTDB,1)]';
% group name translation: specific for rapid tolerance
MWTDB = parseToleranceName(MWTDB);


%% IMPORT
Import_ShaneSpark_Time;
%% CAL: MEAN 
Cal_PlateMean;
%% TABLE OUTPUT: MEAN PER GROUP (N=PLATES)
Cal_GroupMean;
%% MAKE GRAPH
Graph_4dots;
%% MULTI-FACTORIAL ANOVA: STRAIN X PREDOSE X POSTDOSE (N=PLATES)
Stats_mANOVA_NPlate;
save(fullfile(pSave,'data.mat'),'MWTSet');
%% SCATTER PLOT SPEED VS CURVE - BY PLATE
Graph_Scatter_SpeedxCurve;
%% CALCULATE EXP MEAN
save(fullfile(pSave,'data.mat'),'MWTSet');
%%
Cal_ExpMean;
Stats_mANOVA_NExp;
%% GRAPH SCATTER - by exp
Graph_Cond_clusterDotsNExp;
Graph_Cond_clusterDotsNPlate;
%% INDEX -----------------------
Cal_EI;
%% EI GRAPHS
Graph_EI_clusterDotsNPlate;
Graph_EI_clusterDotsNExp;

%% Export
Export_Finish;







