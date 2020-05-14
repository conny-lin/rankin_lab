function CS = CurveStats_v1707r2(DataMeta,MWTDB,msr)
%% get data
data = DataMeta.(msr);
mwtid = DataMeta.mwtid;
%% get group name
a = MWTDB(:,{'mwtid','groupname'});
c = innerjoin(a,DataMeta);
CS.groupnames = c.groupname;
%% get exp name
a = MWTDB(:,{'mwtid','expname'});
c = innerjoin(a,DataMeta);
CS.expnames = c.expname;

%% get descriptpive by plate
S = statsBasicG(data,mwtid,'mwtid'); % make table
a = MWTDB(:,{'mwtid','expname','groupname','strain','rx'});
b = S(:,{'mwtid'});
L = innerjoin(b,a);
CS.dsc_plate = innerjoin(L,S);

%% get descriptive by group
A = CS.dsc_plate;
g = A.groupname;
x = A.mean;
CS.dsc_group = statsBasicG(x,g,'groupname');


%% ANOVA get.dscstat
A = CS.dsc_plate;
x = A.mean;
g = A.groupname;
if numel(x) ~= numel(g); error('check'); end
if numel(unique(g))>1
    [anovatext,anovastats,multstats,T,ST] = anova1_std_v2(x,g);
else
    anovatext = 'no comparison made b/c only 1 group';
    % create desc stats
    T = statsBasic(x,'outputtype','struct');
    % correct naming convention
    T.gnames = unique(g);
    T.SE = T.se;
    T.N = T.n; 
    T.n = [];
end

CS.dscstat = T;
CS.anovatext = anovatext;
CS.anova = ST;
CS.posthoc = multstats;


%% get.statstr(obj)
% descriptive
A = CS.dscstat;
gg = A.group;
A = sortN2first(gg,A);

g = char(strjoinrows(A.group',', '));
n = char(strjoin(num2cellstr(A.n),', '));
m = char(strjoin(num2cellstr(A.mean),', '));
s = char(strjoin(num2cellstr(A.se),', '));
CS.dsctxt = sprintf('G = %s\nN = %s\nmean = %s\nSE = %s',g,n,m,s);
%--------------------------------------------------------------------------



%% get.meanCtrl(obj)
S = CS.dsc_plate;
C = S(ismember(S.rx,'NA'),:);
C.ctrl_name = strjoinrows([C.expname C.strain],'*');
SC = statsBasicG(C.mean,C.ctrl_name,'ctrl_name');
CL = unique(C(:,{'groupname','strain','rx','expname','ctrl_name'}));
SC = innerjoin(CL,SC);
SC.ctrl = SC.mean;
SC.mean = [];
CS.meanCtrl = SC;

%% get.meanExp(obj)
S = CS.dsc_plate;
E = S(~ismember(S.rx,'NA'),:);
E.ctrl_name = strjoinrows([E.expname E.strain],'*');
CS.meanExp = E;

%% Percent ---------------------------------------------------------------


%% get.pctByPlate(obj)
E = CS.meanExp;
C = CS.meanCtrl(:,{'ctrl_name','ctrl'});
% divide against control
A = innerjoin(E,C);
A.pct = (A.mean./A.ctrl)-1;
CS.pct_plate = A;


%% get.pctdscstat(obj)
A = CS.pct_plate;
x = A.pct;
g = A.groupname;
CS.pct_dscstat = statsBasicG(x,g,'groupname');


%% anova of pct stats 
A = CS.pct_plate;
x = A.pct;
g = A.groupname;
if numel(x) ~= numel(g); error('check'); end
if numel(unique(g))>1
    [anovatext,anovastats,multstats,T,ST] = anova1_std_v2(x,g);
else
    anovatext = 'no comparison made b/c only 1 group';
    % create desc stats
    T = statsBasic(x,'outputtype','struct');
    % correct naming convention
    T.gnames = unique(g);
    T.SE = T.se;
    T.N = T.n; 
    T.n = [];
end
CS.pct_dscstat = T;
CS.pct_anova = ST;
CS.pct_posthoc = multstats;



%% get.pct_tstat(obj)
A = CS.pct_plate;
gn = unique(A.groupname);
gn = sortN2first(gn,gn);
T = table;
T.gn = gn;
T.tstat = cell(size(T.gn,1),1);
T.pv = nan(size(T.gn,1),1);
T.pvs = cell(size(T.gn,1),1);
for gi =1:numel(gn)
    a = A.pct(ismember(A.groupname,gn{gi}));
    [~,p,~,stats] = ttest(a,0);
    T.tstat{gi} = sprintf('t(%d) = %.3f',stats.df,stats.tstat);   
    T.pv(gi) = p;
    T.pvs{gi} = print_pvalue(p);
end
CS.pct_ttest_stat = T;

% ttest
a = CS.pct_ttest_stat;
b = strjoinrows([a.gn a.tstat a.pvs],', ');
CS.pct_ttest_str = print_cellstring(b);



% str
CS.statstr = sprintf('*** Descriptive ***\n%s\n\n*** ANOVA ***\n%s\n\n*** ttest 0 ***\n%s\n',...
    CS.dsctxt,CS.anovatext,CS.pct_ttest_str);
end



    
