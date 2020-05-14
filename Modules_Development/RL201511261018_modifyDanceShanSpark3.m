

%%
pHome = cd;
pMWT = D.mwtpath;
% pMWT = PlateInfo.mwtpath;
% pMWT = regexprep(T.pMWT,'/Users/connylin/Dropbox/Lab/MWT_Analysis',pDrive);
addpath('/Users/connylin/Dropbox/RL/Code/Modules/Dance_ShaneSpark3')
MWTSet = Dance_ShaneSpark3(pMWT,'pSave',pHome);

