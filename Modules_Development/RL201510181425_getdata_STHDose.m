%% function RL201510161612_get60mindata

% add function
addpath('/Users/connylin/Dropbox/MATLAB/Function_Library_Public');

%% get path to experiment with dose
% experiments with dose analysis are tagged as "Dose"
pData = '/Volumes/ParahippocampalGyrus/MWT/Data';
[~,~,f,p] = dircontent(pData);
i = regexpcellout(f,'Dose');
pExp = p(i);
expnameT = f(i);
%% check if all has analysis output
% get path to all experiment folders
pA = '/Users/connylin/Dropbox/Lab/MWT_Analysis';
[~,~,f,p] = dircontent(pA);
if sum(ismember(expnameT,f)) ~= numel(expnameT); 
    error('some exp does not have chor analysis files'); 
else
    i = ismember(f,expnameT);
    pExp = p(i);
end

%% create experiment set folder
pDest = '/Users/connylin/Dropbox/Lab/MWT_Analysis_ByExp/Dose_10sISI';
if isdir(pDest) == 0; mkdir(pDest); end


%% create summary
[~,~,fG,pG] = cellfun(@dircontent,pExp,'UniformOutput',0);
pG = celltakeout(pG);
[~,~,fMWT,pMWT] = cellfun(@dircontent,pG,'UniformOutput',0);
% get MWT names and paths
fMWT = celltakeout(fMWT);
pMWT = celltakeout(pMWT);
% make table
T = table;
T.mwtname = fMWT;
[pG,fMWT] = cellfun(@fileparts,pMWT,'UniformOutput',0);
pG = celltakeout(pG);
[pE,fG] = cellfun(@fileparts,pG,'UniformOutput',0);
fG = celltakeout(fG);
pE = celltakeout(pE);
T.groupname = fG;
[~,fE] = cellfun(@fileparts,pE,'UniformOutput',0);
fE = celltakeout(fE);
T.expname = fE;

cd(pDest);
writetable(T,'plate_info.csv');

%% copy to sort into groups
fGU = unique(fG);
for fGUi = 1:numel(fGU)
    % get folder name
    pDestG = [pDest,'/',fGU{fGUi}];
    if isdir(pDestG) == 0; mkdir(pDestG); end
    i = ismember(fG,fGU(fGUi));
    pT = pMWT(i);
    for pTi = 1:numel(pT)
        ps = pT{pTi};
        [~,fn] = fileparts(ps);
        pd = [pDestG,'/',fn];
        if isdir(pd) == 0; mkdir(pd); end
        copyfile(ps,pd,'f')
    end
end



% for x = 2:numel(pExp)
%     ps = pExp{x};
%     [~,fn] = fileparts(ps);
%     pd = [pDest,'/',fn];
%     if isdir(pd) == 0; mkdir(pd); end
%     copyfile(ps,pd,'f')
% end



















