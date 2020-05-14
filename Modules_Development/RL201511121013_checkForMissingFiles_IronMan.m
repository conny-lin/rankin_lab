
%% check 
addpath('/Users/connylin/Dropbox/MATLAB/Function_Library_Public');
mwtnameSearch = '\<\d{8}[_]\d{6}[.]zip\>';
excludeList = {'Analysis'; 'MatlabAnalysis'};
fMWTExcludeList = {'20120330_132929'};

%% get zipped data
pHome = '/Volumes/COBOLT/MWT_Data_Archive_Zip';
[fE,pE,fEf,pEf] = dircontent(pHome);
% get group folder
[~,~,fG,pG] = cellfun(@dircontent,pEf,'UniformOutput',0);
fG = celltakeout(fG);
pG = celltakeout(pG);
[fMWT,pMWT] = cellfun(@dircontent,pG,'UniformOutput',0);
pMWT = celltakeout(pMWT,'multirow');
fMWT = celltakeout(fMWT,'multirow');
% [~,fMWT] = cellfun(@fileparts,pMWT,'UniformOutput',0);
i = regexpcellout(fMWT,mwtnameSearch);
pMWT(~i) = [];
fMWT(~i) = [];
% take out zip from file name
fMWTZip = fMWT;
fMWT = regexprep(fMWT,'[.]zip','');
% [~,fMWT] = cellfun(@fileparts,pMWT,'UniformOutput',0);


%% get MWT data from checking folder
pA = '/Volumes/IRONMAN/MWT_Data_201406';
pMWTA = {};
[~,~,fEA,pEA] = dircontent(pA);
% exclusions
i = regexpcellout(fEA,'catrina');
fEA(i) = []; pEA(i) = [];
[f,p,fG,pG] = cellfun(@dircontent,pEA,'UniformOutput',0);
f = celltakeout(f,'multirow');
p = celltakeout(p,'multirow');
% check if any MWT files underneath folder
i = regexpcellout(f,'(\<\d{8}[_]\d{6}\>)|(\<\d{8}[_]\d{6}[.]zip\>)');
pMWTA = [pMWTA;p(i)];
fG = celltakeout(fG,'multirow');
pG = celltakeout(pG,'multirow');
[f,p,~,~] = cellfun(@dircontent,pG,'UniformOutput',0);
f = celltakeout(f,'multirow');
p = celltakeout(p,'multirow');
i = regexpcellout(f,'(\<\d{8}[_]\d{6}\>)|(\<\d{8}[_]\d{6}[.]zip\>)');
pMWTA = [pMWTA;p(i)];




%% get mwt names
[~,fMWTA] = cellfun(@fileparts,pMWTA,'UniformOutput',0);

% check missing files
i = ~ismember(fMWTA,fMWT);

char(fMWTA(i))

% pMWTA = celltakeout(pMWTA,'multirow');
% fMWTA = celltakeout(fMWTA,'multirow');
% fMWTAU = unique(fMWTA);

%% code stop
return

%%
20140907_163630 2
20120330_132929  
20141224_142717 2
20141224_142717 3
20150310_135917  
20150310_135917  
20150310_150343  
20150310_150343  
20150310_142634  
20150310_142634  
20150310_152112  
20150310_152112  
20150310_144455  
20150310_144455  
20150310_153845  
20150310_153845  
20150310_155724  
20150310_155724  
20150310_144455  
20150310_144455  
20150310_153845  
20150310_155724  
20150310_155724 
%%

i = find(regexpcellout(fMWT,'20120330_132929'))
pMWT{i}


i = find(regexpcellout(fMWTA,'20120330_132929'))
pMWTA{i}


%% move zip MWT files
pHome = '/Volumes/COBOLT/MWT';
pZipA = '/Volumes/COBOLT/MWT_Data_Archive_Zip';
for x = 101:numel(pMWT)
    ps = pMWT{x};
    [pG,pmwt] = fileparts(ps);
    pd = regexprep(pG,pHome,pZipA);
    if isdir(pd) == 0; mkdir(pd); end
    
    movefile(ps,pd,'f')
end

return
%% zip MWT files
pTemp = '/Volumes/COBOLT/MWT_ZipTemp';
reportJump = 100:100:numel(pMWT);
for x = 4606:numel(pMWT)
    if sum(reportJump == x) > 0
        fprintf('Processing %d/%d MWT files\n',x,numel(pMWT));
    end
    pmwt = pMWT{x};
    
    % check if zip already exist
    [pg,fmwt] = fileparts(pmwt);
    a = dircontent(pg,[fmwt,'.zip']);
    
    if isempty(a) == 1
        [mwtcontent,pmwtC] = dircontent(pmwt);
        i = regexpcellout(mwtcontent,'([.]blob)|([.]blobs)|([.]png)|([.]set)|([.]summary)');
        if sum(i) > 1
            ps = pmwtC(i);
            fprintf('zipping [%s]\n',fmwt);
            zip(pmwt,ps);
        end

    else
        fprintf('[%s] already zipped\n',fmwt);
    end
    % delete extra
    [m,pm] = dircontent(pmwt,'*swanlake2all*dat');
    if isempty(m) == 0
        cellfun(@delete,pm);
    end
end




