%% function RL201510161612_get60mindata

% add function
addpath('/Users/connylin/Dropbox/MATLAB/Function_Library_Public');

%% get path to experiment with 60mins recording
pData = '/Volumes/ParahippocampalGyrus/MWT/Data';
[~,~,f,p] = dircontent(pData);
i = regexpcellout(f,'3600s0x0s0s');
pExp = p(i);

%% copy 
pDest = '/Volumes/IRONMAN/MWT_Data_ByExp/AcuteAlcohol_60min';

for x = 2:numel(pExp)
ps = pExp{x};
[~,fn] = fileparts(ps);
pd = [pDest,'/',fn];
if isdir(pd) == 0; mkdir(pd); end
copyfile(ps,pd,'f')
end

%% create plate info
% original experiment name
% plate orientation
% tracker
% expter




















