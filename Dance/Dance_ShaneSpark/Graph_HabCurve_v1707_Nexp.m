function MWTSet = Graph_HabCurve_v1707_Nexp(MWTSet,msrlist)


%% get data
DataTrv = MWTSet.Data_ExpMean; 
MWTDB = MWTSet.MWTDB;
pSave = MWTSet.pSave;
sname = 'Graph_byExp'; % struct save name
% make group graph data
[i,j] = ismember(DataTrv.mwtid, MWTDB.mwtid);
DataTrv.groupname = MWTDB.groupname(j(i));
DataTrv.strain = MWTDB.strain(j(i));
gnu = output_sortN2first(unique(DataTrv.groupname)); % group name
tapu = unique(DataTrv.tap); % get unique tap

G = struct;
header = struct;
header.tap = repmat(tapu,numel(gnu));
for gi = 1:numel(gnu)
    for msri = 1:numel(msrlist)
        i = ismember(DataTrv.groupname,gnu(gi));
        d = table2array(DataTrv(i,msrlist{msri}));
        g = table2array(DataTrv(i,{'tap'}));
        T = grpstatsTable(d,g,'gnameutitle','tap');
  
        % correction (20170606) -- start
        if iscellstr(T.tap) % correct tap to numeric
            T.tap = cellfun(@str2num,T.tap);
        end
        % correction (20170606) -- end

        G.(msrlist{msri}).groupname(gi) = gnu(gi);
        G.(msrlist{msri}).tap(:,gi) = T.tap;
        G.(msrlist{msri}).Mean(:,gi) = T.mean;
        G.(msrlist{msri}).N(:,gi) = T.n;
        G.(msrlist{msri}).SE(:,gi) = T.se;
    end
end
MWTSet.Graph.ByGroupPerExp = G; % 20170821
%--------------------------------------------------------------------------
% Graph 
%--------------------------------------------------------------------------
DataG = MWTSet.Graph.ByGroupPerExp; % 20170821
for msri = 1:numel(msrlist)
    gn = regexprep(DataG.(msrlist{msri}).groupname,'_',' ');
    y = DataG.(msrlist{msri}).Mean;
    x = DataG.(msrlist{msri}).tap;
    e = DataG.(msrlist{msri}).SE;
    msr = msrlist{msri};
    Graph_HabCurveSS(x,y,e,gn,msr,pSave)
end
%--------------------------------------------------------------------------

% Graph excel ouptut %-----------------------------------------------------
statsnames = {'Mean','SE','N'};
GT = {'','','tap'};
D = MWTSet.Graph.ByGroupPerExp.(msrlist{1}).tap(:,1)'; % 20170821
for si  = 1:numel(statsnames)
    sname = statsnames{si};
    A = []; G = {}; M = {};
    for msri = 1:numel(msrlist)
        msr = msrlist{msri};
        A = [A; DataG.(msr).(sname)'];
        G = [G; DataG.(msr).groupname'];
        M = [M; repmat({msr},numel(DataG.(msr).groupname),1)];
    end
    GT = [GT;  [M G repmat({sname},numel(G),1)]];
    D = [D;A];
end
T = table;
T.msr = GT(:,1);
T.groupname = GT(:,2);
T.stats = GT(:,3);
% create header
t = DataG.(msrlist{1}).tap(:,1)';
ts = num2cellstr(t);
a = repmat({'t'},numel(t),1);
colnames = strjoinrows([a ts'],'')';
T1 = array2table(D,'VariableNames',colnames);
T = [T T1];
writetable(T,fullfile(pSave,'Descriptive HabCurve NExp.csv'));

