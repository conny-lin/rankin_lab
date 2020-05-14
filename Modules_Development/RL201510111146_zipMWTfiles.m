
%% check 
addpath('/Users/connylin/Dropbox/MATLAB/Function_Library_Public');
mwtnameSearch = '\<\d{8}[_]\d{6}\>';
pHome = '/Volumes/COBOLT/MWT';
pDeleteArchive = '/Volumes/COBOLT/MWT_dup_to_delete';
excludeList = {'Analysis'; 'MatlabAnalysis'};

%% analysis folder mwt list
% pA = '/Volumes/COBOLT/MWT_dup_to_delete';
% [~,~,fEA,pEA] = dircontent(pA);
% [~,~,fGA,pGA] = cellfun(@dircontent,pEA,'UniformOutput',0);
% pGA = celltakeout(pGA,'multirow');
% [~,~,fMWTA,pMWTA] = cellfun(@dircontent,pGA,'UniformOutput',0);
% pMWTA = celltakeout(pMWTA,'multirow');
% fMWTA = celltakeout(fMWTA,'multirow');
% fMWTAU = unique(fMWTA);

%% get data home
[fE,pE,fEf,pEf] = dircontent(pHome);
% get group folder
[~,~,fG,pG] = cellfun(@dircontent,pEf,'UniformOutput',0);
fG = celltakeout(fG);
pG = celltakeout(pG);
[~,~,~,pMWT] = cellfun(@dircontent,pG,'UniformOutput',0);
pMWT = celltakeout(pMWT,'multirow');
[~,fMWT] = cellfun(@fileparts,pMWT,'UniformOutput',0);
i = regexpcellout(fMWT,mwtnameSearch);
pMWT(~i) = [];
[~,fMWT] = cellfun(@fileparts,pMWT,'UniformOutput',0);

a = tabulate(fMWT);
a(cell2mat(a(:,2))<2,:) = [];


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




