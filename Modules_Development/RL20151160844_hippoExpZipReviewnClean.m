
%% check 
addpath('/Users/connylin/Dropbox/MATLAB/Function_Library_Public');
mwtnameSearch = '(\<\d{8}[_]\d{6}\>)|(\<\d{8}[_]\d{6}[.]zip\>)';
excludeList = {'Analysis'; 'MatlabAnalysis'};
pStd = '/Volumes/COBOLT/MWT';
pCheck = '/Volumes/Public/Data/MWT_20151111/ExpZip';    
pOut = '/Volumes/IRONMAN/unzip';
pA = pOut;

[~,~,~,pEf] = dircontent(pStd);
% get group folder
[~,~,fG,pG] = cellfun(@dircontent,pEf,'UniformOutput',0);
fG = celltakeout(fG);
pG = celltakeout(pG);
[~,~,fMWT,pMWT] = cellfun(@dircontent,pG,'UniformOutput',0);
pMWT = celltakeout(pMWT,'multirow');
fMWT = celltakeout(fMWT,'multirow');
i = regexpcellout(fMWT,'(\<\d{8}[_]\d{6}\>)|(\<\d{8}[_]\d{6}[.]zip\>)');
pMWT(~i) = [];
fMWT(~i) = [];
fMWTD = fMWT;
pMWTD = pMWT;

%% get MWT data from checking folder
[fEAZip,pEAZip] = dircontent(pCheck,'*.zip');
cd(pStd)
i = regexpcellout(fEAZip,'\<[.]\w*');
pEAZip(i) = [];
fEAZip(i) = [];

for x = 1:numel(pEAZip)
    pz = pEAZip{x};
    pzo = regexprep(pz,'[.]zip','');
    if isdir(pOut) == 0; mkdir(pOut); end
    cd(pOut);
    unzip(pz);
    delete(pz)
    cd(pStd)
    tsfNewFiles2SameMWTFolder(pMWTD,pOut)
    rmdir(pOut,'s');
end

%% check consistency of prefix in MWT files
%% check extra MWT files found in D not in archive (or just rezip)

%%

return

%%
i = regexpcellout(fEAZip,'\<13%*');
fEAZip(i)

%%
ps = pEAZip(i);
pd = '/Volumes/IRONMAN/Keep';
fprintf('\nCopying...');
cellfun(@copyfile,ps,cellfunexpr(ps,pd),cellfunexpr(ps,'f'))
fprintf('DONE\n');



%%
[~,e] = cellfun(@fileparts,cellfun(@fileparts,cellfun(@fileparts,a.path,'UniformOutput',0),'UniformOutput',0),'UniformOutput',0)
eU = unique(e)
% char(eU)
t = eU(3:end-1);
val = false(size(fEAZip,1),1);
for x =3:numel(t)
    i = regexpcellout(fEAZip,['\<',t{x},'*']);
    fEAZip(i)
    ps = pEAZip(i);
    pd = '/Volumes/IRONMAN/Keep';
    fprintf('\nCopying...');
    cellfun(@copyfile,ps,cellfunexpr(ps,pd),cellfunexpr(ps,'f'))
    fprintf('Deleting...');
    cellfun(@delete,ps);
    cellfun(@unzip,regexprep(ps,pCheck,pd))
    fprintf('DONE\n');
end

%%
regexprep(ps,pCheck,pd)

%              
%              
% 20130702B_NG_100s30x10s10s_DA609&KP1097
% 20140604C_SJ_100s30x10x10s_goa1        
% 20140616C_SJ_100s30x10s10s             
% 20140616X_SJ                           
% 20140623C_AH_24hrspreexposure          
% 20150719B_JS_24hrPE_DAx609             
% 20150901C_ST_100s30x10s10s

%%
[f,p,~,pG] = cellfun(@dircontent,pEA,'UniformOutput',0);
f = celltakeout(f,'multirow');
p = celltakeout(p,'multirow');
% check if any MWT files underneath folder
i = regexpcellout(f,'(\<\d{8}[_]\d{6}\>)|(\<\d{8}[_]\d{6}[.]zip\>)');
pMWTA = [pMWTA;p(i)];
pG = celltakeout(pG,'multirow');
[f,p,~,~] = cellfun(@dircontent,pG,'UniformOutput',0);
f = celltakeout(f,'multirow');
p = celltakeout(p,'multirow');
i = regexpcellout(f,mwtnameSearch);
pMWTA = [pMWTA;p(i)];
fMWTA = regexprep(cellfun(@fileparts,pMWTA,'UniformOutput',0),'[.]zip','');


%% get mwt names
[~,fMWTA] = cellfun(@fileparts,pMWTA,'UniformOutput',0);

% check missing files
i = ~ismember(fMWTA,fMWT);
if sum(i) ~= numel(i)
fMWTA(i)
disp(char(pMWTA(i)));
fprintf('\n%d/%d files in T not in A\n',sum(i), numel(i));
else
    fprintf('\nAll files accounted for\n');
end

%%
unique(cellfun(@fileparts,cellfun(@fileparts,pMWTA(i),'UniformOutput',0),'UniformOutput',0))
%%


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
