% add function
addpath('/Users/connylin/Dropbox/MATLAB/Function_Library_Public');
pAnalysis = '/Users/connylin/Dropbox/Lab/MWT_Analysis';
pData = '/Volumes/ParahippocampalGyrus/MWT/Data';

%% get path to experiment with dose
% experiments with dose analysis are tagged as "Dose"
[~,~,f,p] = dircontent(pData);
i = regexpcellout(f,'Dose');
pExp = p(i);
expnameT = f(i);
% get MWT dir
[~,~,~,pG] = cellfun(@dircontent,pExp,'UniformOutput',0);
pG = celltakeout(pG);
[~,pMWT] = cellfun(@dircontent,pG,cellfunexpr(pG,'*.zip'),'UniformOutput',0);
pMWT = celltakeout(pMWT);



%% chor Trinity
addpath('/Users/connylin/Dropbox/MATLAB/programs_rankinlab/Library/Modules/Chor');
[pMWTA,fval,chorscript] = chormaster3('Trinity', pMWT,pAnalysis,pData,'old');



%% get path to experiment with dose
% check if all has analysis output
fnameVal = fval{1};
[fn,p] = cellfun(@dircontent,pMWTA,cellfunexpr(pMWTA,fnameVal),'UniformOutput',0);
nFile = cellfun(@numel,fn);
nVal = nFile > 0;
if sum(nVal) ~= numel(pMWTA)
    error('not all MWT has %s output',fnameVal);
else
    fprintf('all %d MWT files have %s outputs\n',numel(pMWTA),fnameVal);
end

% create experiment set folder
pDest = '/Users/connylin/Dropbox/Lab/MWT_Analysis_ByExp/Dose_10sISI';
if isdir(pDest) == 0; mkdir(pDest); end

% find mwt file info
[pG, fMWT] = cellfun(@fileparts,pMWTA,'UniformOutput',0);
[pExp, fG] = cellfun(@fileparts,pG,'UniformOutput',0);
[~, fE] = cellfun(@fileparts,pExp,'UniformOutput',0);

% make table
T = table;
T.mwtname = fMWT;
T.groupname = fG;
T.expname = fE;
cd(pDest);
writetable(T,'plate_info.csv');

% copy to sort into groups
fGU = unique(fG);
for fGUi = 1:numel(fGU)
    % get folder name
    pDestG = [pDest,'/',fGU{fGUi}];
    if isdir(pDestG) == 0; mkdir(pDestG); end
    i = ismember(fG,fGU(fGUi));
    pT = pMWTA(i);
    for pTi = 1:numel(pT)
        ps = pT{pTi};
        [~,fn] = fileparts(ps);
        pd = [pDestG,'/',fn];
        if isdir(pd) == 0; mkdir(pd); end
        copyfile(ps,pd,'f')
    end
end

















