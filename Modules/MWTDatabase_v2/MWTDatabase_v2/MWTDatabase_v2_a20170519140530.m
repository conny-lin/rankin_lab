%% INITIALIZING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc; clear; close all;
addpath('/Users/connylin/Dropbox/Code/Matlab/Library/General');
pM = setup_std(mfilename('fullpath'),'RL','genSave',true);
addpath(pM);
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% SETTINGS     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pNewData = '/Volumes/COBOLT/MWT_New';
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% GET MWT INFO        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

%% group names       +++++++++++++++++++++++++++++++++++++++++++++++++++++++
[groupnames,pG] = cellfun(@dircontent,pExp,'UniformOutput',0);
groupnames = celltakeout(groupnames);
v = cell2mat(regexpi(groupnames,'^[A-Z]{1,}\d{1,}'));
if sum(v) == numel(groupnames) % if names are valid
    fprintf('Validated: group names\n');
else
    fprintf('Error: group names below are incorrect: \n');
    disp(char(groupnames(~v)))
    error('Please correct above exp names');
end
return
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

% MWT names         +++++++++++++++++++++++++++++++++++++++++++++++++++++++
[~,pMWT] = cellfun(@dircontent,celltakeout(pG),'UniformOutput',0);
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

pMWT = celltakeout(pMWT);

return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% %% INITIALIZING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% clc; clear; close all;
% addpath('/Users/connylin/Dropbox/Code/Matlab/Library/General');
% pM = setup_std(mfilename('fullpath'),'RL','genSave',true);
% % addpath(pM);
% 
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
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