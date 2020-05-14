%% NPY
pHome = '/Users/connylin/Dropbox/RL/MWT_Analysis_ByExp/10sISI by Strain';
Ctrl = 'N2';

pData = '/Volumes/COBOLT/MWT';
D = load([pData,'/MWTDatabase.mat']);
MWTDatabase = D.MWTDatabase;
Db = MWTDatabase.mwt;
i = (~cellfun(@isempty,regexp(Db.genotype,'hawaiian')));
p = '/Users/connylin/Dropbox/RL/MWT_Analysis_ByExp/NPY/neuropeptide Y 10sISI/raster';
s = unique([regexprep(dircontent(p),'_400mM','');unique(Db.strain(i,:))]);

i = (ismember(Db.strain,s) | ...
        ismember(Db.strain,'N2')) &...
    ismember(Db.rx,{'NA','400mM'}) &...
    Db.preplate == 100 &...
    Db.tapN == 30 &...
    Db.ISI == 10;
Db(~i,:) = [];
strainList = unique(Db.strain);
strainList(ismember(strainList,'N2')) = [];


for sli = 1:numel(strainList)
    strainT = strainList{sli};
   
    % find target strain
    i = ismember(Db.strain,strainT) &...
        Db.preplate == 100 &...
        Db.tapN == 30 &...
        Db.ISI == 10;
    p = Db.mwtpath(i);
    DbT = Db(i,:);
    % find control from strain exp
    [~,fExp] = cellfun(@fileparts,cellfun(@fileparts,cellfun(@fileparts,p,'UniformOutput',0),'UniformOutput',0),'UniformOutput',0);
    fExpU = unique(fExp);
    i = ismember(Db.expname,fExpU) &...
        ismember(Db.strain,Ctrl);
    DbT = [DbT;Db(i,:)];
    % make sure records are unique
    [~,b] = unique(DbT.mwt_id);
    DbT = DbT(b,:);

    % report on group
    if sum(ismember(unique(DbT.groupname),[strainT,'_400mM']))== 0
        fprintf('\nNo [%s_400mM] -- skip\n',strainT)
    else
        pSave = sprintf('%s/%s',pHome,strainT);
        if isdir(pSave) == 0; mkdir(pSave); end

        pMWT = DbT.mwtpath;
        addpath('/Users/connylin/Dropbox/RL/Code/Modules/Dance_Glee')
        MWTSet = Dance_Glee(pMWT,'pSave',pSave);
    end
end















