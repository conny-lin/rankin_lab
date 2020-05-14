%% function RL201510161612_get60mindata

% add function
addpath('/Users/connylin/Dropbox/MATLAB/Function_Library_Public');

%% get path to experiment with dose
% experiments with dose analysis are tagged as "Dose"
pData = '/Volumes/ParahippocampalGyrus/MWT/Data';
[~,~,f,p] = dircontent(pData);
i = regexpcellout(f,'Dose');
pExp = p(i);
expnameT = f(i);
%% check if all has analysis output
% % get path to all experiment folders
% pA = '/Users/connylin/Dropbox/Lab/MWT_Analysis';
% [~,~,f,p] = dircontent(pA);
% if sum(ismember(expnameT,f)) ~= numel(expnameT); 
%     error('some exp does not have chor analysis files'); 
% else
%     i = ismember(f,expnameT);
%     pExp = p(i);
% end

%% get MWT dir
[~,~,~,pG] = cellfun(@dircontent,pExp,'UniformOutput',0);
pG = celltakeout(pG);
[~,~,~,pMWT] = cellfun(@dircontent,pG,'UniformOutput',0);
pMWT = celltakeout(pMWT);


















