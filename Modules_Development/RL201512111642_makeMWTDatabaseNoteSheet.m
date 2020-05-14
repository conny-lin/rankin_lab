
%% query database
pData = '/Volumes/COBOLT/MWT';
D = load([pData,'/MWTDatabase.mat']);
MWTDB = D.MWTDatabase;
% get targets
Db = MWTDB.mwt;


cd(pData); 
T = table;
T.mwtid = Db.mwt_id;
T.mwtpath = Db.mwtpath;
T.mwtname = Db.mwt;
T.expname = Db.expname;
T.groupname = Db.groupname;
writetable(T,'MWTDB_Note.csv');