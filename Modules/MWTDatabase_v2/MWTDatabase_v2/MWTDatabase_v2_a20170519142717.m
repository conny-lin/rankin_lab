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
[en,pExp] = dircontent(pNewData);
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
[groupnames,pG] = cellfun(@dircontent,pExp,'UniformOutput',0);
groupnames = celltakeout(groupnames);
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
[mwtname,pMWT] = cellfun(@dircontent,celltakeout(pG),'UniformOutput',0);
mwtname = celltakeout(mwtname); % get MWT names
pMWT = celltakeout(pMWT); % get MWT path

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

% update database +++++++++++++++++++++++++++++++++++++++++++++++++++++++++
MWTDB = makeMWTDatabase3(pData,pMWTDB);
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



fprintf('integration of new data completed \n\n');



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% 
% pDataHome = '/Volumes/COBOLT/MWT_New';
% [~,expFolder] = dircontent(pDataHome);
% 
% 
% for expFolderi = 1:numel(expFolder)
%     en = expFolder{expFolderi};
%     [b, c] =dircontent(en);
%     
%     i = regexpcellout(b,'^[.]');
%     groupFolder = c(~i);
%     
%     for groupFolderi = 1:numel(groupFolder)
%        gn = groupFolder{groupFolderi};
%        
%         
%         [b, c] =dircontent(gn);
%         
%         i = regexpcellout(c,'[.]zip$');
%         zipfilename = c(i);
%         
%         for zi = 1:numel(zipfilename)
%             n = zipfilename{zi};
%             fprintf('unzip: %s',n)
% 
%             unzip(n);
%             
%             delete(n);
%             fprintf(' (done)\n')
%             
%         end
%     end
%     
%     
% end