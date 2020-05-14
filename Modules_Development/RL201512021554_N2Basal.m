%% remake MWT database
pData = '/Volumes/COBOLT/MWT';
pAH = '/Users/connylin/Dropbox/RL/MWT_Analysis_ByExp';
addpath('/Users/connylin/Dropbox/MATLAB/Function_Library_Public');
pSave = sprintf('%s/Basal 5mins',pAH);
if isdir(pSave) == 0; mkdir(pSave); end

% load database
D = load([pData,'/MWTDatabase.mat']);
MWTDatabase = D.MWTDatabase;
% determine
rcT = {'0s2x0s1320s'};
Db = MWTDatabase.mwt;
% assay groups for rc
unique(Db.groupname(ismember(Db.rc,rcT)))



%% list of other rc for basal activities
% '0s2x0s1320s'
% '300s0x0s0s' - age effect not done
% '3600s0x0s0s' - done
% '600s0x0s0s'
% '900s0x0s0s

%% create condition id index, manually add description
a = unique(Db.groupname(ismember(Db.rc,rcT)));
b = regexpcellout(a,'_','split');
s = b(:,1);
cond = b(:,2:end);
% find unique condition
cond = reshape(cond,numel(cond),1);
cond(cellfun(@isempty,cond)) = [];
condU = unique(cond);
t = table;
t.cond_id = (1:numel(condU))';
t.cond = condU;
cd(pHome)
writetable(t,'condition_id.csv');


%% group list for 300s

grouplist = {...
    {'N2 4d 300s' {'N2' 'N2_400mM'}};...
    {'BZ142 4d 300s' {'N2' 'N2_400mM' 'BZ142' 'BZ142_400mM'}};...
    {'DA609 4d 300s' {'N2' 'N2_400mM' 'DA609' 'DA609_400mM'}};...
    {'N2 3d 24h exposure' {'N2_E3d24h0mM_R1h_T4d0mM' 'N2_E3d24h0mM_R1h_T4d200mM' 'N2_E3d24h200mM_R1h_T4d0mM' 'N2_E3d24h200mM_R1h_T4d200mM'}};...
    {'DA609 3d 24h exposure' {'DA609_E3d24h0mM_R1h_T4d0mM' 'DA609_E3d24h0mM_R1h_T4d200mM' 'DA609_E3d24h200mM_R1h_T4d0mM' 'DA609_E3d24h200mM_R1h_T4d200mM'}};...
    {'MT14390' {'MT14390_E3d24h0mM_R1h_T4d0mM' 'MT14390_E3d24h0mM_R1h_T4d200mM' 'MT14390_E3d24h200mM_R1h_T4d0mM' 'MT14390_E3d24h200mM_R1h_T4d200mM'}};...
    {'MT14480' {'MT14480_E3d24h0mM_R1h_T4d0mM' 'MT14480_E3d24h0mM_R1h_T4d200mM' 'MT14480_E3d24h200mM_R1h_T4d0mM' 'MT14480_E3d24h200mM_R1h_T4d200mM'}};...
    {'RB1025' {'RB1025_E3d24h0mM_R1h_T4d0mM' 'RB1025_E3d24h0mM_R1h_T4d200mM' 'RB1025_E3d24h200mM_R1h_T4d0mM' 'RB1025_E3d24h200mM_R1h_T4d200mM'}};...
};
% 'N2' 'N2_400mM'
% 'N2' 'N2_400mM' 'BZ142' 'BZ142_400mM'
% 'N2' 'N2_400mM' 'DA609' 'DA609_400mM'
% 'N2' 'N2_400mM' 'NM1968' 'NM1968_400mM'
% 'N2_E3d24h0mM_R1h_T4d0mM' 'N2_E3d24h0mM_R1h_T4d200mM' 'N2_E3d24h200mM_R1h_T4d0mM' 'N2_E3d24h200mM_R1h_T4d200mM'

% 'N2_E3d24h0mM_R1h_E4d0mM_T5d0mM'
% 'N2_E3d24h200mM_R1h_E4d0mM_T5d0mM'
% 'N2_E3d24h0mM_R1h_T4d0mM' 'N2_E3d24h0mM_R1h_T4d200mM' 'N2_E3d24h200mM_R1h_T4d0mM' 'N2_E3d24h200mM_R1h_T4d200mM' 'DA609_E3d24h0mM_R1h_T4d0mM' 'DA609_E3d24h0mM_R1h_T4d200mM' 'DA609_E3d24h200mM_R1h_T4d0mM' 'DA609_E3d24h200mM_R1h_T4d200mM'
% 'MT14390_E3d24h0mM_R1h_T4d0mM' 'MT14390_E3d24h0mM_R1h_T4d200mM' 'MT14390_E3d24h200mM_R1h_T4d0mM' 'MT14390_E3d24h200mM_R1h_T4d200mM'
% 'MT14480_E3d24h0mM_R1h_T4d0mM' 'MT14480_E3d24h0mM_R1h_T4d200mM' 'MT14480_E3d24h200mM_R1h_T4d0mM' 'MT14480_E3d24h200mM_R1h_T4d200mM'
% 'RB1025_E3d24h0mM_R1h_T4d0mM' 'RB1025_E3d24h0mM_R1h_T4d200mM' 'RB1025_E3d24h200mM_R1h_T4d0mM' 'RB1025_E3d24h200mM_R1h_T4d200mM'


%% group list archive
% group list for 3600s0x0s0s ----------------------------------------------
% grouplist = {...
%     {'N2 dose 200 400' {'N2' 'N2_200mM' 'N2_400mM'}};...
%     {'N2 age 0mM' {'N2' 'N2_3d' 'N2_5d'}};...
%     {'N2 age 400mM' {'N2_400mM_3d' 'N2_400mM' 'N2_400mM_5d'}};...
%     {'BZ142' {'N2' 'N2_400mM' 'BZ142' 'BZ142_400mM'}};...
%     {'KC565' {'N2' 'N2_400mM' 'KC565' 'KC565_400mM'}};...
%     {'MT14390' {'N2' 'N2_400mM' 'MT14390' 'MT14390_400mM'}};...
%     {'MT14480' {'N2' 'N2_400mM' 'MT14480' 'MT14480_400mM'}};...
%     {'NM1630' {'N2' 'N2_400mM' 'NM1630' 'NM1630_400mM'}};...
%     {'NM1968' {'N2' 'N2_400mM' 'NM1968' 'NM1968_400mM'}};...
%     {'RB1025 200mM' {'N2' 'N2_200mM' 'RB1025' 'RB1025_200mM'}};...
%     {'RB1025' {'N2' 'N2_400mM' 'RB1025' 'RB1025_400mM'}};...
%     {'VC1038' {'N2' 'N2_400mM' 'VC1038' 'VC1038_400mM'}};...
%     {'VC922' {'N2' 'N2_400mM' 'VC922' 'VC992_400mM'}}...
% };



%% all the groups for  3600s0x0s0s
% 'N2'
% 'N2_200mM'
% 'N2_3d'
% 'N2_400mM'
% 'N2_400mM_3d'
% 'N2_400mM_5d'
% 'N2_5d'
% 'BZ142'
% 'BZ142_400mM'
% 'KC565'
% 'KC565_400mM'
% 'MT14390'
% 'MT14390_400mM'
% 'MT14480'
% 'MT14480_400mM'
% 'NM1630'
% 'NM1630_400mM'
% 'NM1968'
% 'NM1968_400mM'
% 'RB1025'
% 'RB1025_200mM'
% 'RB1025_400mM'
% 'VC1038'
% 'VC1038_400mM'
% 'VC922'
% 'VC992_400mM'



%% filter
addpath('/Users/connylin/Dropbox/RL/Code/Modules/Dance_DrunkMoves');
% get targets
Db = MWTDatabase.mwt;
% find exp with N2_400mM group
for gli = 1:numel(grouplist)
    i = ismember(Db.rc,rcT) & ...
        ismember(Db.groupname,grouplist{gli}{1,2});
    pMWT = Db.mwtpath(i);
    setname = grouplist{gli}{1,1};
    pSaveA = sprintf('%s/%s',pSave,setname);
    if isdir(pSaveA) == 0; mkdir(pSaveA); end
    MWTSet = Dance_DrunkMoves(pMWT,'pSave',pSaveA);
    
end

return














