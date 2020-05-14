% function RL201510021354_takeFreqFigureOut(p)
% based on RL201509251335_postanalysis_STH


% add function for shared folder
addpath('/Users/connylin/Dropbox/MATLAB/Function_Library_Public');

% get experiment info in the analysis folder
pData = '/Users/connylin/Dropbox/Lab/MWT_Analysis';
[~,~,En,pE] = dircontent(pData);
iExpVal = true(size(En));
%% get after liquid transfer
% get exp date
a = regexpcellout(En,'^\d{8}','match');
date = cellfun(@str2num,a);
datelimit = 20120127;
i = date <= datelimit;
iExpVal(i) = false;

% get 10sISI 
term = '100s30x10s';
i = regexpcellout(En,term);
iExpVal(~i) = false;

%% get only alcohol worm group experiments
% get expname
expter = regexpcellout(En,'(?<=^\d{8}[A-Z]_)[A-Z]{2}','match');
% unique(expter)
i = ismember(expter,{'PO'});
iExpVal(i) = false;

%% export
EnT = En(iExpVal);
disp(EnT);

t = cell2table(EnT,'VariableNames',{'expname'});
cd('/Users/connylin/Dropbox/Lab/MWT_ExpSet');
writetable(t,'alcohol_10sISI.csv');

%% run Dance only on N2
pET = pE(iExpVal);

%% get group
[~,~,gn,pg] = cellfun(@dircontent,pET,'UniformOutput',0);
gn = celltakeout(gn);
pg = celltakeout(pg);
pg(~ismember(gn,{'N2','N2_400mM'})) = []; % remove not N2 files
addpath('/Users/connylin/Dropbox/MATLAB/Programs_RankinLab/Library/Dance_ShaneSpark3');
MWTSet = Dance_ShaneSpark3(pg,2);   




