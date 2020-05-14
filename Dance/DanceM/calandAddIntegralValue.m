function Tbl = calandAddIntegralValue(X1,G1,Y1,Tbl,si)

% calculate under the curve for each plate
nr = max(unique(X1));
mwtidu = unique(G1);
nc = max(mwtidu);
B = nan(nr,nc);
i = sub2ind([nr,nc],X1,G1);
B(i) = Y1;
X = repmat([1:30]',1,nr);
A = diff(X);
y = B(1:end-1,:);
y2 = diff(B)./2;
y = y+y2;
Y = nansum(y);
D1 = table;
D1.mwtid = [1:max(mwtidu)]';
D1.area = Y';

D2 = CurveStats;
D2.curve = D1.area;
D2.mwtid = D1.mwtid;
D2.MWTDB = M;
[anovatxt,T] = anova(D2);
gn = regexprep(T.gnames,'_400mM','');

a = sortN2first(gn,T.N)';
a = num2cellstr(a);
a = strjoin(a,', ');
Tbl.N{si} = a;

b = (sortN2first(gn,T.mean).*100)-100;
Tbl.wt(si) = b(1);
Tbl.mut(si) = b(2);

d = b(2) - b(1);
Tbl.d_mut(si) =d;

a = regexp(anovatxt,'(?<=p (<|=) )\d{1,}[.]\d{1,}','match');
a = cellfun(@str2num,a);
Tbl.p(si) = a;

A = D2.pctEtohByPlate;
x1 = A.pct_ctrl(ismember(A.groupname,'N2_400mM'));
x2 = A.pct_ctrl(ismember(A.groupname,[strain,'_400mM']));
e = effectsize_cohend(x1,x2);
Tbl.ES(si) = e;   