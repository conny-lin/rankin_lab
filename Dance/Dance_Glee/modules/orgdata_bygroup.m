function A = orgdata_bygroup(D)
% D = MWTSet.Data.trv;

A = struct;
msru = fieldnames(D);
for mi = 1:numel(msru)
    msr = msru{mi};
    d = D.(msr);
    p = d.pMWT;
    a = parseMWTinfo(p);
    gnu = unique(a.groupname);
    for gi = 1:numel(gnu)
        gn = gnu{gi};
        i = ismember(a.groupname,gn);

        fnu = fieldnames(d);
        c = struct;
        for fi = 1:numel(fnu)
            fn = fnu{fi};
            c.(fn) = d.(fn)(:,i);
        end
        A.(msr).(gn) = c;
    end
end
