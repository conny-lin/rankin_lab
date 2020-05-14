%% INITIALIZING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc; clear; close all;
addpath('/Users/connylin/Dropbox/Code/Matlab/Library/General');
pM = setup_std(mfilename('fullpath'),'RL','genSave',true);
addpath(pM);
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% SETTINGS     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pNewData = '/Volumes/COBOLT/MWT_New';
pDataBase = '/Volumes/COBOLT/MWT';
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% PROCESS NEW FILES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% experiment names ++++++++++++++++++++++++++++++++++++++++++++++++++++++++
[~,pExp] = dircontent(pNewData); % survey for files
pExp = takeout_tempfiles(pExp); % take out temp files
en = cellfun(@fileparts,pExp); % get exp name
v = regexpcellout(en,'^\d{8}[A-Z][_][A-Z]{2}[_]\d{1,}s\d{1,}x\d{1,}s\d{1,}s'); % name validator
if sum(v) == numel(en) % if names are valid
    fprintf('Validated: experiment names\n');
else
    fprintf('Error: experiment names below are incorrect: \n');
    disp(char(en(~v)))
    error('Please correct above exp names');
end
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

% group names       +++++++++++++++++++++++++++++++++++++++++++++++++++++++
[~,pG] = cellfun(@dircontent,pExp,'UniformOutput',0);
pG = celltakeout(pG); % take out group names paths
pG = takeout_tempfiles(pG); % take out temp files
groupnames = cellfun(@fileparts,pG); % get file name
v = cell2mat(regexpi(groupnames,'^[A-Z]{1,}\d{1,}'));
if sum(v) == numel(groupnames) % if names are valid
    fprintf('Validated: group names\n');
else
    fprintf('Error: group names below are incorrect: \n');
    disp(char(groupnames(~v)))
    error('Please correct above group names');
end

%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

% MWT names         +++++++++++++++++++++++++++++++++++++++++++++++++++++++
[~,pMWT] = cellfun(@dircontent,celltakeout(pG),'UniformOutput',0);
pMWT = celltakeout(pMWT); % get MWT path
pMWT = takeout_tempfiles(pMWT); % take out temp files
mwtname = cellfun(@fileparts,pMWT); % get file name
% check if zipped
z = regexpcellout(mwtname,'[.]zip$');
pMWTzipped = pMWT(z);
fprintf('%d/%d MWT files zipped\n',sum(z),numel(pMWT));
if sum(z) > 0
   fprintf('unzipping:\n');
    for zi = 1:numel(pMWTzipped)
        processIntervalReporter(numel(pMWTzipped),1,'unzipping',zi)
        p = pMWTzipped{zi};
        p2 = fileparts(p(1:end-4));
        unzip(p,p2);
        delete(p);
    end
    fprintf('unzip; done\n');
end
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%% INTEGRATE TO MAIN DATABASE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% move new files to database ++++++++++++++++++++++++++++++++++++++++++++++
[en,pExp] = dircontent(pNewData);
for i = 1:numel(pExp)
    ps = pExp{i};
    pd = pDataBase;
    movefile(ps,pd);
end
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

%% update database ++++++++++++++++++++++++++++++++++++++++++++++++++++++++
MWTDB = makeMWTDatabase3(pDataBase,pDataBase);
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



fprintf('integration of new data completed \n\n');



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

