function MWTDB = DanceM_MWTDBUpdate(pMWT,pDataBase)
% update MWTDB from database

Db = load(pDataBase); Db = Db.MWTDB.text;
i = ismember(Db.mwtpath,pMWT);
if sum(i) ~= numel(pMWT)
   error('database missing input pMWT info');
end
[i,j] = ismember(pMWT,Db.mwtpath);
MWTDB = Db(j(i),:); clear Db;
MWTDB.mwtid_db = MWTDB.mwtid;
MWTDB.mwtid = [1:size(MWTDB,1)]';