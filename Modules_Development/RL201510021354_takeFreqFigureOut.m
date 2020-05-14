function RL201510021354_takeFreqFigureOut(p)
% based on RL201509251335_postanalysis_STH

%% move all freq figures into one file
pStrainFiles = [p,'/Graph Habituation curves by strain'];
[~,~,~,pG] = dircontent(pStrainFiles);
[~,pF] = cellfun(@dircontent,pG,'UniformOutput',0);
pF = celltakeout(pF);
pF(cellfun(@isempty,regexp(pF,'RevFreq'))) = [];

% define save folder
pSave = [fileparts(pStrainFiles),'/Graph Hab curves Freq'];
if isdir(pSave) == 0; mkdir(pSave); end
%% create designation paths
for x = 1:numel(pF)
pS = pF{x};
pSf = fileparts(pS);
pD = regexprep(pS,pSf,pSave);
copyfile(pS,pD)
end



