function pMiss = checkForMissingFiles(pStd,pCheck)
%% checkForMissingFiles(pStd,pCheck)
% pStd = not zipped files
% pCheck = zip files

%% check 
% addpath('/Users/connylin/Dropbox/MATLAB/Function_Library_Public');
mwtnameSearch = '(\<\d{8}[_]\d{6}\>)|(\<\d{8}[_]\d{6}[.]zip\>)';
% excludeList = {'Analysis'; 'MatlabAnalysis'};
% pHome = '/Volumes/IRONMAN/MWT_Data_Zip';
% pStd = '/Volumes/COBOLT/MWT';
% pCheck = '/Volumes/IRONMAN/MWT_Data_Zip';

%% get zipped data
fprintf('getting MWT folder data from first path\n');
[~,~,~,pEf] = dircontent(pStd);
% get group folder
[~,~,~,pG] = cellfun(@dircontent,pEf,'UniformOutput',0);
pG = celltakeout(pG);
[fMWT,pMWT] = cellfun(@dircontent,pG,'UniformOutput',0);
pMWT = celltakeout(pMWT,'multirow');
fMWT = celltakeout(fMWT,'multirow');
pMWT(~regexpcellout(fMWT,mwtnameSearch)) = [];
% take out zip from file name
[~,fMWT] = cellfun(@fileparts,pMWT,'UniformOutput',0);


%% get MWT data from checking folder
fprintf('getting MWT zip data from second path\n');
pMWTA = {};
[~,~,~,pEA] = dircontent(pCheck);
[f,p,~,pG] = cellfun(@dircontent,pEA,'UniformOutput',0);
f = celltakeout(f,'multirow');
p = celltakeout(p,'multirow');
% check if any MWT files underneath folder
pMWTA = [pMWTA;p(regexpcellout(f,mwtnameSearch))];
pG = celltakeout(pG,'multirow');
[f,p,~,~] = cellfun(@dircontent,pG,'UniformOutput',0);
p = celltakeout(p,'multirow');
pMWTA = [pMWTA;p(regexpcellout(celltakeout(f,'multirow'),mwtnameSearch))];
% fMWTA = regexprep(cellfun(@fileparts,pMWTA,'UniformOutput',0),'[.]zip','');


%% get mwt names
[~,fMWTA] = cellfun(@fileparts,pMWTA,'UniformOutput',0);
% check missing files
i = ~ismember(fMWTA,fMWT);
if sum(i) ~= numel(i)
    char(fMWTA(i))
    disp(char(pMWTA(i)));
    fprintf('\n%d/%d files in T not in A\n',sum(i), numel(i));
else
    fprintf('\nAll files accounted for\n');
end
disp(char(unique(cellfun(@fileparts,cellfun(@fileparts,pMWTA(i),'UniformOutput',0),'UniformOutput',0))))

pMiss = pMWTA(i);

return

% pMWTA = celltakeout(pMWTA,'multirow');
% fMWTA = celltakeout(fMWTA,'multirow');
% fMWTAU = unique(fMWTA);
% 
% %%
% i = find(regexpcellout(fMWT,'20150910_174545'))
% pMWT{i}
% 
% 
% %% move zip MWT files
% pHome = '/Volumes/COBOLT/MWT';
% pZipA = '/Volumes/COBOLT/MWT_Data_Archive_Zip';
% for x = 101:numel(pMWT)
%     ps = pMWT{x};
%     [pG,pmwt] = fileparts(ps);
%     pd = regexprep(pG,pHome,pZipA);
%     if isdir(pd) == 0; mkdir(pd); end
%     
%     movefile(ps,pd,'f')
% end
% 
% 
% %% zip MWT files
% pTemp = '/Volumes/COBOLT/MWT_ZipTemp';
% reportJump = 100:100:numel(pMWT);
% for x = 4606:numel(pMWT)
%     if sum(reportJump == x) > 0
%         fprintf('Processing %d/%d MWT files\n',x,numel(pMWT));
%     end
%     pmwt = pMWT{x};
%     
%     % check if zip already exist
%     [pg,fmwt] = fileparts(pmwt);
%     a = dircontent(pg,[fmwt,'.zip']);
%     
%     if isempty(a) == 1
%         [mwtcontent,pmwtC] = dircontent(pmwt);
%         i = regexpcellout(mwtcontent,'([.]blob)|([.]blobs)|([.]png)|([.]set)|([.]summary)');
%         if sum(i) > 1
%             ps = pmwtC(i);
%             fprintf('zipping [%s]\n',fmwt);
%             zip(pmwt,ps);
%         end
% 
%     else
%         fprintf('[%s] already zipped\n',fmwt);
%     end
%     % delete extra
%     [m,pm] = dircontent(pmwt,'*swanlake2all*dat');
%     if isempty(m) == 0
%         cellfun(@delete,pm);
%     end
% end
% 
% 
% 
% 
