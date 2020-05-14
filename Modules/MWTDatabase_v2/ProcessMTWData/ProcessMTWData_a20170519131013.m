%% INITIALIZING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc; clear; close all;
addpath('/Users/connylin/Dropbox/Code/Matlab/Library/General');
pM = setup_std(mfilename('fullpath'),'RL','genSave',true);
% addpath(pM);
return
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


pDataHome = '/Volumes/COBOLT/MWT_New';
[~,expFolder] = dircontent(pDataHome);


for expFolderi = 1:numel(expFolder)
    en = expFolder{expFolderi};
    [b, c] =dircontent(en);
    
    i = regexpcellout(b,'^[.]');
    groupFolder = c(~i);
    
    for groupFolderi = 1:numel(groupFolder)
       gn = groupFolder{groupFolderi};
       
        
        [b, c] =dircontent(gn);
        
        i = regexpcellout(c,'[.]zip$');
        zipfilename = c(i);
        
        for zi = 1:numel(zipfilename)
            n = zipfilename{zi};
            fprintf('unzip: %s',n)

            unzip(n);
            
            delete(n);
            fprintf(' (done)\n')
            
        end
    end
    
    
end