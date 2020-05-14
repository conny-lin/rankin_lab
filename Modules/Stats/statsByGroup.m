function G = statsByGroup(Data,msrlist,timevar,grpname)

% % get group info ++++
% M = parseMWTinfo(MWTDB.mwtpath(Data.mwtid));
% a = table;
% a.groupname = M.groupname;
% Data = [a Data];
% % ------------------


gnu = unique(Data.(grpname));
G = struct;

for msri = 1:numel(msrlist)
    
    msr = msrlist{msri};
        
    G.(msr).(grpname) = gnu;
    
    for gi = 1:numel(gnu)
        
        gn = gnu{gi};
        i = ismember(Data.(grpname),gn);       
        
        x = Data.(msr)(i);
        g = Data.(timevar)(i);
        a = statsBasicG(x,g);
        
        if gi==1
            X = nan(size(a,1),numel(gnu));
            Y = X;
            E = X;
            N = X;
        end
        
        X(:,gi) = a.gname;
        N(:,gi) = a.n;
        Y(:,gi) = a.mean;
        E(:,gi) = a.se;
        
    end
    
    G.(msr).X = X;
    G.(msr).N = N;
    G.(msr).Y = Y;
    G.(msr).E = E;
end