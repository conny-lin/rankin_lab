%% 60sISI

%% user input variables
pHome = '/Users/connylin/Dropbox/RL/MWT_Analysis_ByExp/60sISI by Strain';
pData = '/Volumes/COBOLT/MWT';
ISI = 60;
tapN = 30;
Ctrl = 'N2';

%% process variables
if isdir(pHome) == 0; mkdir(pHome); end
D = load([pData,'/MWTDatabase.mat']);
MWTDatabase = D.MWTDatabase;

%% get strains
Db = MWTDatabase.mwt;
i = ...
    ismember(Db.rx,{'400mM'}) &...
    Db.preplate == 100 &...
    Db.tapN == tapN &...
    Db.ISI == ISI;
s = unique(Db.groupname(i));
s = regexprep(s,'[_]400mM','');
s(ismember(s,'N2')) = [];
% make sure all strains have corresponding no alcohol group
i = ismember(s,Db.groupname);
if sum(i) ~= numel(i)
    error('some strains do not have 0mM control group');
end
a = cellfun(@strcat,s,cellfunexpr(s,'_400mM'),'UniformOutput',0);
grouplist = [s;a];
strainlist = s;


%% make sure each strains have N2 control, if not, manually enter control name
ctrlNameList = cell(size(strainlist));
for x = 1:numel(strainlist)
    strain = strainlist{x};
    i = ismember(Db.groupname,[strain,'_400mM']);
    explist = Db.expname(i);
    i = ismember(Db.expname,explist);
    if sum(ismember(Db.groupname(i),{'N2_400mM'})) == 0 ||...
        sum(ismember(Db.groupname(i),{'N2'})) == 0
        warning('%s has no N2 and N2_400mM',strain)    
    else
        ctrlNameList{x} = 'N2';
    end    
end
if sum(cellfun(@isempty,ctrlNameList)) > 0
   error('enter control name manually'); 
%    ctrlNameList(ismember(strainlist,'VG202')) = {'NM1968'};
%    ctrlNameList(ismember(strainlist,'VG302')) = {'BZ142'};
end


%% run Glee
addpath('/Users/connylin/Dropbox/RL/Code/Modules/Dance_Glee')
for si = 1:numel(strainlist)
    strainT = strainlist{si};
    refStrain = ctrlNameList{si};
    grouplist = {refStrain;[refStrain,'_400mM'];strainT;[strainT,'_400mM']};
    explist = unique(Db.expname(ismember(Db.groupname,{strainT;[strainT,'_400mM']})));
    i = ismember(Db.expname,explist) &...
        ismember(Db.groupname,grouplist) &...
        Db.preplate == 100 &...
        Db.tapN == tapN &...
        Db.ISI == ISI;
    pMWT = Db.mwtpath(i);
    pSave = sprintf('%s/%s',pHome,strainT);
    if isdir(pSave) == 0; mkdir(pSave); end
    MWTSet = Dance_Glee(pMWT,'pSave',pSave,'refStrain',refStrain);
end


return





%% SUMMARY OF RESULTS

%% produce summary for alcohol effect
R = cell(size(strainlist));
for si = 1:numel(strainlist)
    strain = strainlist{si};
    refStrain = ctrlNameList{si};
    str = sprintf('%s/%s/Dance_Glee/AlcoholEffectSummary.csv',pHome,strain);
    R{si} = readtable(str);
end

% reconstruct table
msrlist = {'RevDur','RevFreq','RevSpeed'};
assaylist = {'Initial','HabLevel','HabRate_integral'};
A = table; nRow = 1;
for mi = 1:numel(msrlist)
    msrname = msrlist{mi};
    for ai = 1:numel(assaylist)
        assayname = assaylist{ai};
        for si = 1:numel(strainlist)
            r = R{si};
            i = ismember(r.msr,msrname) & ismember(r.assays,assayname);
            a = table;
            a.msr = {msrname};
            a.assays = {assayname};
            fn = fieldnames(r);
            refstrain = fn{3};
            teststrain = fn{4};
            a.refstrain = {refstrain};
            a.teststrain = {teststrain};
            a.refalcoholeffect = r.(refstrain)(i);
            a.testalcoholeffect = r.(teststrain)(i);
            a.ref_test_STHeffect = r.test_vs_ref_0mM(i);
            a.test_alcoholeffect = r.ref_alcoholeffect_vs_ref(i);
            A = [A;a];
        end
    end
end
A.test_alcoholeffect(cellfun(@isempty,A.test_alcoholeffect)) = {'same trend'};
writetable(A,sprintf('%s/AlcoholEffectStrainSummary.csv',pHome));



%% summarize
pSave = [pHome,'/strain sig summary']; if isdir(pSave) == 0; mkdir(pSave); end

msrlist = {'RevFreq','RevDur','RevSpeed'};
assaylist = {'Initial','HabLevel','HabRate_integral'};
refEffectlist = unique(A.refalcoholeffect);
effectnamelist = unique(A.test_alcoholeffect);
effectnamelist(ismember(effectnamelist,{'--'})) = [];
refstrain = 'N2';


fid = fopen([pSave,'/summary.txt'],'w');
for mi = 1:numel(msrlist)
for tei = 1:numel(effectnamelist)
for ai = 1:numel(assaylist)
for rei = 1:numel(refEffectlist)
effectname = effectnamelist{tei};
msrname = msrlist{mi};
assayname = assaylist{ai};
refEffect = refEffectlist{rei};

i = ismember(A.test_alcoholeffect,effectname) & ...
    ismember(A.msr,msrname) &...
    ismember(A.assays,assayname) &...
    ismember(A.refalcoholeffect,refEffect) & ...
    ismember(A.refstrain,refstrain);
if sum(i) > 0
    A(i,:);

    s = A.teststrain(i);
    [i,j] = ismember(s,Db.strain);
    
    T = table;
    T.strain = Db.strain(j);
    T.genotype = Db.genotype(j);
    str = sprintf('%s %s %s %s',msrname,assayname,refEffect,effectname);
    writetable(T,sprintf('%s/%s.txt',pSave,str));
    
    fprintf(fid,'%s:\n',str);
    for x = 1:size(T,1)
       fprintf(fid,'%s, %s\n',T.strain{x},T.genotype{x}); 
    end
    fprintf(fid,'\n');
end
end
end
end
end
fclose(fid);


%% expset
fid = fopen([pSave,'/expset.txt'],'w');
expset = MWTDatabase.expset;
for x = 1:size(expset,1)
    fprintf(fid,'%s:\n',expset.expset{x});
    a = expset.genes{x};
    for xx = 1:numel(a)
        fprintf(fid,'%s\n',a{xx});
    end
    fprintf(fid,'\n');
end
fclose(fid);


%% summarize - table
pSave = [pHome,'/strain sig summary']; if isdir(pSave) == 0; mkdir(pSave); end

msrlist = {'RevFreq','RevDur','RevSpeed'};
assaylist = {'Initial','HabLevel','HabRate_integral'};
refEffectlist = unique(A.refalcoholeffect);
effectnamelist = unique(A.test_alcoholeffect);
effectnamelist(ismember(effectnamelist,{'--'})) = [];
refstrain = 'N2';
%% create strainlist table by group
expset = MWTDatabase.expset;
strainU = unique(A.teststrain);
[i,j] = ismember(strainU,Db.strain);
T = table;
T.strain = strainU;
T.genotype = Db.genotype(j);
GenotypeList = T;

%%
expset = MWTDatabase.expset;
T = table;
for x = 1:size(expset,1)
    a = expset.genes{x};
    if size(a,1) == 1
        a = a';
    end
    e = repmat(expset.expset(x),size(a,1),1);
    t = table;
    t.set = e;
    t.gene = a;
    if size(T,1) == 0;
        T = t;
    else
        T = [T;t];
    end
end
SetList = T;

%% create strain and genotype list
T = table;
T.Set = cell(size(GenotypeList,1),1);
T = [T GenotypeList];
for x = 1:size(SetList,1)
    str = sprintf('^%s*',SetList.gene{x});
    i = regexpcellout(T.genotype,str);
    T.Set(i) = SetList.set(x);
    
end
% manually enter loose ends
T.Set(ismember(T.strain,'DA650')) = {'neuropeptide Y'};
T.Set(ismember(T.strain,'JPS383')) = {'BK'};
T.Set(ismember(T.strain,'RM3389')) = {'neuroligin'};
T = sortrows(T,{'Set'});
Legend = T;
Legend.teststrain = Legend.strain;
Legend.strain = [];

%%
D = innerjoin(Legend,A);
D.genotype = regexprep(D.genotype,',',';');
writetable(D,sprintf('%s/Summary with gene groups.csv',pSave));

%% separate by set
setlist = unique(D.Set);
for si = 1:numel(setlist)
   setname = setlist{si};
   d = D(ismember(D.Set,setname),:);
   writetable(d,sprintf('%s/Summary with gene groups - %s.csv',pSave,setname));

end

%%
% fid = fopen([pSave,'/summary.txt'],'w');
for mi = 1:numel(msrlist)
for ai = 1:numel(assaylist)
for tei = 1:numel(effectnamelist)
for rei = 1:numel(refEffectlist)
effectname = effectnamelist{tei};
msrname = msrlist{mi};
assayname = assaylist{ai};
refEffect = refEffectlist{rei};

i = ismember(A.test_alcoholeffect,effectname) & ...
    ismember(A.msr,msrname) &...
    ismember(A.assays,assayname) &...
    ismember(A.refalcoholeffect,refEffect) & ...
    ismember(A.refstrain,refstrain);
if sum(i) > 0
    A(i,:)
    
    error('stop')
    s = A.teststrain(i);
    [i,j] = ismember(s,Db.strain);
    
    T = table;
    T.strain = Db.strain(j);
    T.genotype = Db.genotype(j);
    str = sprintf('%s %s %s %s',msrname,assayname,refEffect,effectname);
    writetable(T,sprintf('%s/%s.txt',pSave,str));
    
%     fprintf(fid,'%s:\n',str);
%     for x = 1:size(T,1)
%        fprintf(fid,'%s, %s\n',T.strain{x},T.genotype{x}); 
%     end
%     fprintf(fid,'\n');
end
end
end
end
end
% fclose(fid);



return













%%
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















