%% remake MWT database
pData = '/Volumes/COBOLT/MWT';
pAH = '/Users/connylin/Dropbox/RL/MWT_Analysis_ByExp';
addpath('/Users/connylin/Dropbox/MATLAB/Function_Library_Public');
pSave = sprintf('%s/Rapid Tolerance/Summary',pAH);
if isdir(pSave) == 0; mkdir(pSave); end

% load database
D = load([pData,'/MWTDatabase.mat']);
MWTDatabase = D.MWTDatabase;
% determine
rcT = {'300s0x0s0s'};
Db = MWTDatabase.mwt;
% assay groups for rc
fprintf('groups found in rc[%s]:\n',char(rcT));
disp(unique(Db.groupname(ismember(Db.rc,rcT))))


%% group list for 300s

grouplist = {...
    {'N2' {'N2_E3d24h0mM_R1h_T4d0mM' 'N2_E3d24h0mM_R1h_T4d200mM' 'N2_E3d24h200mM_R1h_T4d0mM' 'N2_E3d24h200mM_R1h_T4d200mM'}};...
    {'DA609' {'DA609_E3d24h0mM_R1h_T4d0mM' 'DA609_E3d24h0mM_R1h_T4d200mM' 'DA609_E3d24h200mM_R1h_T4d0mM' 'DA609_E3d24h200mM_R1h_T4d200mM'}};...
    {'MT14390' {'MT14390_E3d24h0mM_R1h_T4d0mM' 'MT14390_E3d24h0mM_R1h_T4d200mM' 'MT14390_E3d24h200mM_R1h_T4d0mM' 'MT14390_E3d24h200mM_R1h_T4d200mM'}};...
    {'MT14480' {'MT14480_E3d24h0mM_R1h_T4d0mM' 'MT14480_E3d24h0mM_R1h_T4d200mM' 'MT14480_E3d24h200mM_R1h_T4d0mM' 'MT14480_E3d24h200mM_R1h_T4d200mM'}};...
    {'RB1025' {'RB1025_E3d24h0mM_R1h_T4d0mM' 'RB1025_E3d24h0mM_R1h_T4d200mM' 'RB1025_E3d24h200mM_R1h_T4d0mM' 'RB1025_E3d24h200mM_R1h_T4d200mM'}};...
};

%% filter
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
%     MWTSet = Dance_DrunkMoves_RapidTolerance(pMWT,'pSave',pSaveA);
    load(sprintf('%s/Dance_DrunkMoves_RapidTolerance/Dance_DrunkMoves_RapidTolerance.mat',pSaveA))
    
    pS = [pSaveA,'/worm clump picture']; if isdir(pS) == 0; mkdir(pS); end
    D = MWTSet.Data_Plate.goodnumber;
    A = MWTSet.Info.VarIndex;

    gu = unique(D.groupname);
    for gi = 1:numel(gu)
        % find plate with max number of worms
        d = D(D.groupname == gu(gi),:);
        [~,i] = max(d.mean);
        % get source path
        p = A.pmwt{d.pmwt(i)};
        [~,a] = dircontent(p,'*.png');
        ps = a{1};
        % make des path
        % get mwt name 
        mwtn = A.mwtname{d.mwtname(i)};
        % get group name
        gn = A.groupname{d.groupname(i)};
        % translate group name
        gn = regexprep(gn,'(E3d24h)|(R1h_)|(T4d)','');
        % create ps
        pd = sprintf('%s/Max N %s %s.png',pS,gn,mwtn);
        copyfile(ps,pd);
    end
end



%% list of other rc for basal activities
% '0s2x0s1320s'
% '300s0x0s0s' - age effect not done
% '3600s0x0s0s' - done
% '600s0x0s0s'
% '900s0x0s0s

%% group list for 300s
% 
% grouplist = {...
%     {'N2' {'N2_E3d24h0mM_R1h_T4d0mM' 'N2_E3d24h0mM_R1h_T4d200mM' 'N2_E3d24h200mM_R1h_T4d0mM' 'N2_E3d24h200mM_R1h_T4d200mM'}};...
%     {'DA609' {'DA609_E3d24h0mM_R1h_T4d0mM' 'DA609_E3d24h0mM_R1h_T4d200mM' 'DA609_E3d24h200mM_R1h_T4d0mM' 'DA609_E3d24h200mM_R1h_T4d200mM'}};...
%     {'MT14390' {'MT14390_E3d24h0mM_R1h_T4d0mM' 'MT14390_E3d24h0mM_R1h_T4d200mM' 'MT14390_E3d24h200mM_R1h_T4d0mM' 'MT14390_E3d24h200mM_R1h_T4d200mM'}};...
%     {'MT14480' {'MT14480_E3d24h0mM_R1h_T4d0mM' 'MT14480_E3d24h0mM_R1h_T4d200mM' 'MT14480_E3d24h200mM_R1h_T4d0mM' 'MT14480_E3d24h200mM_R1h_T4d200mM'}};...
%     {'RB1025' {'RB1025_E3d24h0mM_R1h_T4d0mM' 'RB1025_E3d24h0mM_R1h_T4d200mM' 'RB1025_E3d24h200mM_R1h_T4d0mM' 'RB1025_E3d24h200mM_R1h_T4d200mM'}};...
% };
% 
% 
% % 'N2' 'N2_400mM'
% % 'N2' 'N2_400mM' 'BZ142' 'BZ142_400mM'
% % 'N2' 'N2_400mM' 'DA609' 'DA609_400mM'
% % 'N2' 'N2_400mM' 'NM1968' 'NM1968_400mM'
% % 'N2_E3d24h0mM_R1h_T4d0mM' 'N2_E3d24h0mM_R1h_T4d200mM' 'N2_E3d24h200mM_R1h_T4d0mM' 'N2_E3d24h200mM_R1h_T4d200mM'
% 
% % 'N2_E3d24h0mM_R1h_E4d0mM_T5d0mM'
% % 'N2_E3d24h200mM_R1h_E4d0mM_T5d0mM'
% % 'N2_E3d24h0mM_R1h_T4d0mM' 'N2_E3d24h0mM_R1h_T4d200mM' 'N2_E3d24h200mM_R1h_T4d0mM' 'N2_E3d24h200mM_R1h_T4d200mM' 'DA609_E3d24h0mM_R1h_T4d0mM' 'DA609_E3d24h0mM_R1h_T4d200mM' 'DA609_E3d24h200mM_R1h_T4d0mM' 'DA609_E3d24h200mM_R1h_T4d200mM'
% % 'MT14390_E3d24h0mM_R1h_T4d0mM' 'MT14390_E3d24h0mM_R1h_T4d200mM' 'MT14390_E3d24h200mM_R1h_T4d0mM' 'MT14390_E3d24h200mM_R1h_T4d200mM'
% % 'MT14480_E3d24h0mM_R1h_T4d0mM' 'MT14480_E3d24h0mM_R1h_T4d200mM' 'MT14480_E3d24h200mM_R1h_T4d0mM' 'MT14480_E3d24h200mM_R1h_T4d200mM'
% % 'RB1025_E3d24h0mM_R1h_T4d0mM' 'RB1025_E3d24h0mM_R1h_T4d200mM' 'RB1025_E3d24h200mM_R1h_T4d0mM' 'RB1025_E3d24h200mM_R1h_T4d200mM'

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










