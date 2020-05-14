function export_trv_raw(pSave,D)
% D is grouped MWT trv data
% save raw trv output
pSaveA = [pSave,'/Data trv Raw'];
if isdir(pSaveA) == 0; mkdir(pSaveA); end

msru = fieldnames(D);
for mi = 1:numel(msru)
    msr = msru{mi};
    fnu = fieldnames(D.(msr));
    fnu(ismember(fnu,{'pMWT','mwtname'})) = [];   
    DB = parseMWTinfo(D.(msr).pMWT);
    T = table;
    T.expname = DB.expname;
    T.groupname = DB.groupname;
    T.mwtname = DB.mwtname;
    for fi = 1:numel(fnu)
        fn = fnu{fi};
        s = D.(msr).(fn)';
        s = array2table(s);
        s = [T s];
        writetable(s,sprintf('%s/%s %s.csv',pSaveA,msr,fn));
    end
end