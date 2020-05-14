function MWTSet = Dance_RapidTolerance_v1707(pMWT,pSave,varargin)
%% Dance_DrunkMoves_RapidTolerance
% updated: 201603270930


%% DEFAULTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('\n** Running %s **\n',mfilename);
% pDataBase = '/Volumes/COBOLT/MWT/MWTDB.mat';
timeStartSet = 240;
timeEndSet = 300;
msrlist= {'speed','curve'};
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% varargin %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[InputStruct] = vararginProcessor_v1707(whos,varargin)
MWTSet.InputStruct = InputStruct;
return
vararginProcessor;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%% output
[MWTSet,pSave] = DanceM_MWTSetStd_v1707(pMWT,mfilename('fullpath'),varargin,pSave);


%% PREPARE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
MWTSet = DanceM_chor_v1707(MWTSet,'DrunkPostureOnly');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% group name translation: specific for rapid tolerance


%% IMPORT
[MWTSet,Data] = DanceM_importchordata(MWTSet,...
                'chorname','drunkposture2',...
                'timerange',[timeStartSet timeEndSet]);

return
%% CAL: MEAN 
% Cal_PlateMean;
MWTSet = DanceM_datatransform(MWTSet);

MWTSet = DanceM_exportcsv(MWTSet);

DanceM_Stats(MWTSet)

DanceM_Graphs(MWTSet)


%% OLD
%
Cal_GroupMean; % TABLE OUTPUT: MEAN PER GROUP (N=PLATES)
Graph_4dots; % MAKE GRAPH
Stats_mANOVA_NPlate; % MULTI-FACTORIAL ANOVA: STRAIN X PREDOSE X POSTDOSE (N=PLATES)
save(fullfile(pSave,'data.mat'),'MWTSet');
Graph_Scatter_SpeedxCurve; %% SCATTER PLOT SPEED VS CURVE - BY PLATE
save(fullfile(pSave,'data.mat'),'MWTSet'); %% CALCULATE EXP MEAN
Cal_ExpMean;
Stats_mANOVA_NExp;
% %% GRAPH SCATTER - by exp
% Graph_Cond_clusterDotsNExp;
% Graph_Cond_clusterDotsNPlate;
% EI -----------------------
Cal_EI;
% %% EI GRAPHS
% Graph_EI_clusterDotsNPlate;
% Graph_EI_clusterDotsNExp;

% Export
Export_Finish;







