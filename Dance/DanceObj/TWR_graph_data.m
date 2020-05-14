function R = TWR_graph_data(D)

y = D.Mean;
e = D.SE;
x = D.tap;
gn = D.groupname;
D1(:,:,1) = x;
D1(:,:,2) = y;
D1(:,:,3) = e;
% sort by N2
gns = sortN2first(gn',gn');
[i,j] = ismember(gn,gns);
D1 = D1(:,j,:);
X = D1(:,:,1);
Y = D1(:,:,2);
E = D1(:,:,3);
% store
R = struct;
R.X = X;
R.Y = Y;
R.E = E;
R.gn = gn;
