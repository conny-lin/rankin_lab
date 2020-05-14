function [pMWTzip,MWTfn,pG] = getMWTzipped(pData)

%% find available files
p = getalldir(pData);
[~,b] = cellfun(@dircontent,p,cellfunexpr(p,'*.zip'),'UniformOutput',0);
b(cellfun(@isempty,b)) = [];
b = celltakeout(b);
[c,d] = cellfun(@fileparts,b,'UniformOutput',0);
pMWTzip = b;
MWTfn = d;
pG = c;