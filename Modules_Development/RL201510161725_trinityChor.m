%% function RL201510161725_trinityChor
pAnalysis = '/Volumes/IRONMAN/MWT_Analysis_ByExp';
pData = '/Volumes/IRONMAN/MWT_Data_ByExp';
addpath('/Users/connylin/Dropbox/MATLAB/Function_Library_Public');
%% pMWTS
pG = {'/Volumes/IRONMAN/MWT_Data_ByExp/AcuteAlcohol_60min/N2';
        '/Volumes/IRONMAN/MWT_Data_ByExp/AcuteAlcohol_60min/N2_400mM'};
%%
[~,p] = cellfun(@dircontent,pG,'UniformOutput',0);
pMWTS = celltakeout(p);
addpath('/Users/connylin/Dropbox/MATLAB/programs_rankinlab/Library/Modules/Chor');
pMWTA = chormaster3('Trinity', pMWTS,pAnalysis,pData,'old');

