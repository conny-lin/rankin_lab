function changeGroupName(oldname,newname)
%% groupname change
% load alternate names
% cd(fileparts(pSave))
% T = readtable('groupname_alternate.csv');
% gntarget = {'T4dxh0mM_T5d0mM' 'T4dxh400mM_T5d0mM' '400mM_400mM'};
% T(~ismember(T.groupname,gntarget),:) = [];


%% defaults
pData = '/Volumes/COBOLT/MWT';
pZipArchive = '/Volumes/IRONMAN/RL_MWT_Data_Zip';

%% process input
T = table;
T.groupname = oldname;
T.groupname_alternative = newname;

%% find paths to files for name change
% load database
D = load([pData,'/MWTDB.mat']);
MWTDB = D.MWTDB.text;
i = ismember(MWTDB.groupname,T.groupname);
MWTDB(~i,:) = [];
pMWT = MWTDB.mwtpath;

%% make sure same names are found in zip archive
a = regexprep(pMWT,pData,pZipArchive);
for mwti = 1:numel(pMWT)
    p = a{mwti};
    [pf,fn] = fileparts(p);
    f = dircontent(pf,[fn,'.zip']);
    if isempty(f)
        disp(p)
        error('no zip file')
    end
end
pMWTzip = a;

%% combine strain name with new rx name
[i,j] = ismember(MWTDB.rx,T.groupname);
groupname_new = T.groupname_alternate(j(i));


%% make new paths
pMWTnew = cell(size(pMWT));
for mwti = 1:numel(pMWT)
   pmwt = pMWT{mwti};
   pMWTnew{mwti} = regexprep(pmwt,MWTDB.groupname{mwti},groupname_new{mwti});
end
a = parseMWTinfo(pMWTnew);

%% check before proceed
b = [MWTDB.groupname a.groupname];
disp(b);
display 'changing name like above')
c = input('are you sure to preceed? (y=1 n=0)');
if c~=0
    display 'abort';
    return
end

%% change group folder name in database
pGold = cellfun(@fileparts,pMWT,'UniformOutput',0);
pGnew = cellfun(@fileparts,pMWTnew,'UniformOutput',0);
for gi = 1:numel(pGold)
    ps = pGold{gi};
    pd = pGnew{gi};
    if isdir(ps)==1
        movefile(ps,pd,'f')
    end

end

%% change group folder name in archive
pGold = cellfun(@fileparts,pMWTzip,'UniformOutput',0);
pGnew = cellfun(@fileparts,pMWTnew,'UniformOutput',0);
pGnew = regexprep(pGnew,pData,pZipArchive);
for gi = 1:numel(pGold)
    ps = pGold{gi};
    pd = pGnew{gi};
    if isdir(ps)==1
        movefile(ps,pd,'f')
    end
end