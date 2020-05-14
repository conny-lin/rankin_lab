function A = orgdata_bygroup_byexp(D)
% D = MWTSet.Data.trv;

A = struct;
msru = fieldnames(D);
for mi = 1:numel(msru)
    msr = msru{mi};
    d = D.(msr);
    fnu = fieldnames(d);

    p = d.pMWT;
    DB = parseMWTinfo(p);

    gnu = unique(DB.groupname);
    enu = unique(DB.expname);
    for gi = 1:numel(gnu)
        for ei = 1:numel(enu)
            gn = gnu{gi};
            en = enu{ei};
            a = char(regexp(en,'\<\d{8}[A-Z]{1}[_][A-Z]{2}[_]\d+s\d+x\d+s\d+s','match'));
            enName = ['E',a];
            i = ismember(DB.groupname,gn) & ismember(DB.expname,en);
            if sum(i) > 1
                c = struct;
                for fi = 1:numel(fnu)
                    fn = fnu{fi};
                    c.(fn) = d.(fn)(:,i);
                end
                A.(msr).(gn).(enName) = c;
            end
        end
    end
end
