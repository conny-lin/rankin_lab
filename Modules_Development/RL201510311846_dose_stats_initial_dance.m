%% HYPOTHESIS: 
% reversal duration: all doses shortens reversal duration
% response latency: higher doses (400mM+) increases response latency
% reversal probability: probability started to decrease from 300mM < use
% bethoven
% separate data to:
%     1-pre tap: use the same analysis as basal activity (not done)
%     after tap responses: changed within 0.4s (two frames) 
%     2-after tap: + response latency
%     3-post tap response: same as basal activity


%% paths & settings
% exp info path
pHome = '/Users/connylin/Dropbox/Lab/MWT_Analysis_ByExp/Dose_10sISI';
pSave = [pHome,'/Stats'];
% function path
addpath('/Users/connylin/Dropbox/MATLAB/Function_Library_Public');
addpath('/Users/connylin/Dropbox/MATLAB/programs_rankinlab/Library/Modules/rasterPlot_colorSpeed');
addpath('/Users/connylin/Dropbox/MATLAB/programs_rankinlab/Library/Modules/Dance_ShaneSpark3');
pAnalysis = '/Users/connylin/Dropbox/Lab/MWT_Analysis';



%% load global data
plateInfo = readtable(sprintf('%s/plate_info.csv',pHome));



% generate full list of pMWTA
pMWTA = cell(size(plateInfo,1),1);
for x = 1:size(plateInfo,1)
    pMWTA{x} = sprintf('%s/%s/%s/%s',pAnalysis,plateInfo.expname{x}, plateInfo.groupname{x}, plateInfo.mwtname{x});
end


MWTSet = Dance_ShaneSpark3(pMWTA,3)













return















