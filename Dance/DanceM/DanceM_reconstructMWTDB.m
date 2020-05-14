function [MWTSet,MWTDb] = DanceM_reconstructMWTDB(MWTSet)

pMWT = MWTSet.PATHS.pMWT;
MWTDb = parseMWTinfo(pMWT);
mwtname = MWTSet.Info.VarIndex.mwtname;
[i,j] = ismember(mwtname,MWTDb.mwtname);
if sum(i) ~= numel(mwtname)
   error('input pMWT not enough for result mwtname'); 
else
    MWTDb = MWTDb(j(i),:);
end
% reconstruct numbering
if ~isequal(MWTDb.mwtname,mwtname); 
    error('bad matching'); 
else
    MWTDb.mwtid = (1:size(MWTDb,1))';
end
MWTSet.Info.MWTDB = MWTDb;