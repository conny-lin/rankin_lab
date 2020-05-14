function strainNames = DanceM_load_strainInfo(strainlist)

%% load strain genotype list
pDataBase = '/Volumes/COBOLT/MWT';
load(fullfile(pDataBase,'MWTDB.mat')); % load MWTDB
strainNames = MWTDB.strain;

% if no strain list given
if nargin==0
    p = '/Users/connylin/Dropbox/RL/RL Pub PhD Dissertation/Chapters/3-Genes/3-Results/0-Data/10sIS by strains';
    strainlist = dircontent(p);
end

% get strain info
% load('/Users/connylin/Dropbox/Code/Matlab/Library RL/Modules/MWTDatabase/StrainNames.mat');
strainNames(~ismember(strainNames.strain,strainlist),:) = [];

if nargin==0
   strainNames = sortrows(strainNames,{'strain'}); 
end


