%% summarize the robustness of alcohol effct on N2 10sISI

%% user input variables
pHome = '/Users/connylin/Dropbox/RL/MWT_Analysis_ByExp/N2 10sISI effect by exp graph/RevFreq_byeffect';
prefix = 'RevFreq N2 ';
[~,pCat] = dircontent(pHome);
T = table;
for pci = 1:numel(pCat)
    p = pCat{pci};
    [~,fCat] = fileparts(p);
    a = dircontent(p);
    a = regexprep(a,'([.]eps)','');
    a = regexprep(a,prefix,'');
    t = table;
    t.cat = repmat(cellstr(fCat),numel(a),1);
    t.expdate = cellfun(@str2num,regexpcellout(a,'\<\d{8}(?=[A-Z][_])','match'));
    t.tracker = regexpcellout(a,'(?<=\<\d{8})[A-Z]','match');
    t.expter = regexpcellout(a,'(?<=\<\d{8}[A-Z][_])[A-Z]{2}','match');
    T = [T;t];
end
cd(pHome);
writetable(T,'Effect analysis.csv');


return

%%
ISI = 10;
tapN = 30;
Ctrl = 'N2';
xstart = 83;

%% process variables
if isdir(pHome) == 0; mkdir(pHome); end
D = load([pData,'/MWTDatabase.mat']);
MWTDatabase = D.MWTDatabase;

%% get targets
Db = MWTDatabase.mwt;
% find exp with N2_400mM group
i = ...
    ismember(Db.groupname,{'N2_400mM'}) &...
    Db.preplate == 100 &...
    Db.tapN == tapN &...
    Db.ISI == ISI;
expT = unique(Db.expname(i));
% get exp with both N2 and N2_400mM group
i = ...
    ismember(Db.expname,expT) &...
    ismember(Db.groupname,{'N2'});
expT2 = unique(Db.expname(i));
% get paths to target MWT files
i = ...
    ismember(Db.expname,expT2) &...
    ismember(Db.groupname,{'N2','N2_400mM'});
% get target
Db(~i,:) = [];
% get unique experiments list
expU = unique(Db.expname);
% validate targets
% numel(expU)
% unique(DbT.groupname)

%% run Glee
addpath('/Users/connylin/Dropbox/RL/Code/Modules/Dance_Glee')
grouplist = {'N2','N2_400mM'};
for ei = xstart:numel(expU)
    fprintf('\n\n**processing %d/%d exp\n',ei,numel(expU));
    expname = expU{ei};
%     
%     refStrain = ctrlNameList{ei};
%     grouplist = {refStrain;[refStrain,'_400mM'];strainT;[strainT,'_400mM']};
%     explist = unique(Db.expname(ismember(Db.groupname,{strainT;[strainT,'_400mM']})));
    i = ismember(Db.expname,expname) &...
        ismember(Db.groupname,grouplist) &...
        Db.preplate == 100 &...
        Db.tapN == tapN &...
        Db.ISI == ISI;
    pMWT = Db.mwtpath(i);
    pSave = sprintf('%s/%s',pHome,expname);
    if isdir(pSave) == 0; mkdir(pSave); end
    MWTSet = Dance_Glee_Showmance(pMWT,'pSave',pSave);
    if isempty(MWTSet.Data.Raw.MWTfn) ~= 1
        % transfer graph files
        ps = sprintf('%s/Dance_Glee_Showmance/Graph HabCurve/RevFreq N2.eps',pSave);
        pd = sprintf('%s/RevFreq/RevFreq N2 %s.eps',pGraphTsf,expname);
        if isdir(fileparts(pd))== 0; mkdir(fileparts(pd)); end
        copyfile(ps,pd);
        ps = sprintf('%s/Dance_Glee_Showmance/Graph HabCurve/RevDur N2.eps',pSave);
        pd = sprintf('%s/RevDur/RevDur N2 %s.eps',pGraphTsf,expname);
        if isdir(fileparts(pd))== 0; mkdir(fileparts(pd)); end
        copyfile(ps,pd);
        ps = sprintf('%s/Dance_Glee_Showmance/Graph HabCurve/RevSpeed N2.eps',pSave);
        pd = sprintf('%s/RevSpeed/RevSpeed N2 %s.eps',pGraphTsf,expname);
        if isdir(fileparts(pd))== 0; mkdir(fileparts(pd)); end
        copyfile(ps,pd);
    end
end


return



