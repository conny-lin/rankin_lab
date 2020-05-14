function RL201510211328_trinityChor_function(pMWT)

pAnalysis = '/Users/connylin/Dropbox/Lab/MWT_Analysis';
pData = '/Volumes/ParahippocampalGyrus/MWT/Data';
addpath('/Users/connylin/Dropbox/MATLAB/Function_Library_Public');
%% pMWTS
% pG = {'/Volumes/IRONMAN/MWT_Data_ByExp/AcuteAlcohol_60min/N2';
%         '/Volumes/IRONMAN/MWT_Data_ByExp/AcuteAlcohol_60min/N2_400mM'};
%%
% [~,p] = cellfun(@dircontent,pG,'UniformOutput',0);
% pMWTS = celltakeout(p);
addpath('/Users/connylin/Dropbox/MATLAB/programs_rankinlab/Library/Modules/Chor');

pMWTA = chormaster3('Trinity', pMWT,pAnalysis,pData,'old');

end
