%% RL201509251328_Analyze_STH10sISI_Alcohol
addpath('/Users/connylin/Dropbox/MATLAB/Function_Library_Public');

%% Get paths for target experiemnts
p1 = '/Users/connylin/Dropbox/Lab/MWT_Analysis';
[~,~,~,pExp] = dircontent(p1);
[~,~,Gn,pG] = cellfun(@dircontent,pExp,'UniformOutput',0);
Gn = celltakeout(Gn);
pG = celltakeout(pG);
if numel(Gn) ~= numel(pG); error('conversion error'); end
pE = cellfun(@fileparts,pG,'UniformOutput',0);
[~,En] = cellfun(@fileparts,pE,'UniformOutput',0);

%% find target group names
val = true(size(Gn));

val(cellfun(@isempty,regexp(Gn,'400mM'))) = false;
val(~cellfun(@isempty,regexp(Gn,'N2'))) = false;
val(~cellfun(@isempty,regexp(Gn,'N2_400mM'))) = false;
val(~cellfun(@isempty,regexp(Gn,'(NoLid)|(NoFood)|(Test)'))) = false;

a = cell2mat(cellfun(@size,regexp(Gn,'_'),'UniformOutput',0));
val(a(:,2) ~= 1) = false;

val(cellfun(@isempty,regexp(En,'100s30x10s10s'))) = false;
val(~cellfun(@isempty,regexp(En,'Training'))) = false;
val(~cellfun(@isempty,regexp(En,'training'))) = false;


strain_400mM = unique(Gn(val));

a = regexpcellout(strain_400mM,'_','split');
a = unique(a(:,1));
groupNameTarget = sortrows([a;strain_400mM]);

% find path to target groupnames
pGT = pG(val);

% get experiment name
pET = cellfun(@fileparts,pGT,'UniformOutput',0);
[~,ETn] = cellfun(@fileparts,pET,'UniformOutput',0);
ETnU = unique(ETn);

pETU = unique(pET);


%% create group paths for Dance
[~,~,Gn,pG] = cellfun(@dircontent,pETU,'UniformOutput',0);
Gn = celltakeout(Gn);
pG = celltakeout(pG);
if numel(Gn) ~= numel(pG); error('conversion error'); end
pE = cellfun(@fileparts,pG,'UniformOutput',0);
[~,En] = cellfun(@fileparts,pE,'UniformOutput',0);

val = true(size(Gn));
% unique(Gn)
val(~cellfun(@isempty,regexp(Gn,'(NoLid)|(NoFood)|(Test)'))) = false;
unique(Gn(val))
pGT = pG(val);


%% run Dance
addpath('/Users/connylin/Dropbox/MATLAB/Programs_RankinLab/Library/Dance_ShaneSpark3');
Dance_ShaneSpark3(pGT,2)