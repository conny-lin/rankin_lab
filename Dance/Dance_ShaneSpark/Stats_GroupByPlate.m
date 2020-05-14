function G = Stats_GroupByPlate(DataTrv,MWTDB,msrlist)

if nargin<=2
msrlist = {'RevFreq','RevSpeed','RevDur'};
end
% make group graph data
% DataTrv = MWTSet.Raw;
% MWTDB = MWTSet.MWTDB;
[i,j] = ismember(DataTrv.mwtid,MWTDB.mwtid);
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

    G.(msrlist{msri}).groupname(gi) = gnu(gi);
    G.(msrlist{msri}).tap(:,gi) = T.tap;
    G.(msrlist{msri}).Mean(:,gi) = T.mean;
    G.(msrlist{msri}).N(:,gi) = T.n;
    G.(msrlist{msri}).SE(:,gi) = T.se;
end
end
