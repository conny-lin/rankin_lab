function MWTDB = MWTDatabase_v2(pNewData,pDataBase)
%% INSTRUCTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% data must be prepared according to this instruction:
% https://www.evernote.com/l/ADe-7p_cJCJPmae-azp-MFU988pObC454iQ
% file format must be like this: 20130409A_CL_100s30x10s10s_slo1
% 
% pNewData = path of folder that contains new data.
% pDataBase = path of folder that contain the MWT database.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% INITIALIZING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% clc; clear; close all;
addpath('/Users/connylin/Dropbox/Code/Matlab/Library/General');
% pM = setup_std(mfilename('fullpath'),'RL','genSave',true);
% addpath(pM);
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% SETTINGS     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin==0
    pNewData = '/Volumes/COBOLT/MWT_New';
    pDataBase = '/Volumes/COBOLT/MWT';
end
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Clean database tempfiles
[~,p] = dircontent(pDataBase);
takeout_tempfiles(p); % take out temp files


%% PROCESS NEW FILES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% experiment names ++++++++++++++++++++++++++++++++++++++++++++++++++++++++
[~,pExp] = dircontent(pNewData); % survey for files
if isempty(pExp)
    fprintf('no new data found\n');
   return 
end
[pExp,en] = takeout_tempfiles(pExp); % take out temp files
v = regexpcellout(en,'^\d{8}[A-Z][_][A-Z]{2}[_]\d{1,}s\d{1,}x\d{1,}s\d{1,}s'); % name validator
if sum(v) == numel(en) % if names are valid
    fprintf('Validated: experiment names\n');
else
    fprintf('Error: experiment names below are incorrect: \n');
    disp(char(en(~v)))
    error('Please correct above exp names');
end

%% check if experiment already exist in database
[fnD,~] = dircontent(pDataBase);
i = ismember(en,fnD);
if any(i)
    warning('some experiments already exists in database, skip those:');
    disp(en(i)); % display exp
    pExp(i) = []; % exclude exp
    if isempty(pExp)
        disp('no new experiments to add');
        return
    end
end
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

%% group names       +++++++++++++++++++++++++++++++++++++++++++++++++++++++
[~,pG] = cellfun(@dircontent,pExp,'UniformOutput',0);
pG = celltakeout(pG); % take out group names paths
[pG,groupnames] = takeout_tempfiles(pG); % take out temp files
a = regexpi(groupnames,'^[A-Z]{1,}\d{1,}'); % validate names
a(cellfun(@isempty,a)) = {0};
v = cell2mat(a);
numel(v)
if sum(v) == numel(groupnames) % if names are valid
    fprintf('Validated: group names\n');
else
    fprintf('Error: group names below are incorrect: \n');
    disp(char(groupnames(~v))); % display
    error('Please correct above group names');
end

%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

% MWT names         +++++++++++++++++++++++++++++++++++++++++++++++++++++++
[~,pMWT] = cellfun(@dircontent,pG,'UniformOutput',0);
pMWT = celltakeout(pMWT); % get MWT path
[pMWT,~] = takeout_tempfiles(pMWT); % take out temp files
mwtname = cellfun(@dircontent,pG,'UniformOutput',0);
mwtname = celltakeout(mwtname); % get MWT path
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
    movefile(ps,pd,'f');
end
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

%% update database ++++++++++++++++++++++++++++++++++++++++++++++++++++++++
MWTDB = makeMWTDatabase3(pDataBase,pDataBase);
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

fprintf('integration of new data completed \n\n');



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
