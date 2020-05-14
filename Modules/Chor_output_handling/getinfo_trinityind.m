function A = getinfo_trinityind(MWTDB)


%% get information from trinity files
pMWT = MWTDB.mwtpath;

a = cellfun(@dircontent,pMWT,cellfunexpr(pMWT,'trinity_worm*'),'UniformOutput',0);
% create pMWT legend
n = cellfun(@numel,a);
A = cell(size(a,1),1);
for ni = 1:numel(n)
   A{ni} = repmat(pMWT(ni),n(ni),1);
end
A = celltakeout(A);
% A = parseMWTinfo(A);

%%
B = table;
B.mwtpath = A;
C = innerjoin(B,MWTDB(:,{'mwtid','mwtpath','groupname'}));
A = C;
% parse wormid, start and end time
a = celltakeout(a);
a1 = regexpcellout(a,'\d{1,}','match');
A.wormid = cellfun(@str2num,a1(:,1));
% enter in table
A.t1 = cellfun(@str2num,a1(:,2));
A.t2 = cellfun(@str2num,a1(:,3));
A.filename = a;
