function IS_stats(DataMeta,MWTDB,pSave)
gns = MWTDB.groupname(DataMeta.mwtid);

msrlist = {'curve','speedbm','speed'};
for msri = 1:numel(msrlist)
    msr = msrlist{msri};
    [strain,rx] = parse_groupname(gns);
    rx(cellfun(@isempty,rx)) = {'0mM'};
    G = [{strain} {rx}];
    gvar = {'strain','dose'};
    Y = DataMeta.(msr);
    [s,T] = anovan_std(Y,G,gvar,pSave,'suffix',msr);
    writetable(T,fullfile(pSave,sprintf('%s.csv', msr)));
end
