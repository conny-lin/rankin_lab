
%% group name translation: specific for rapid tolerance
MWTDB = MWTSet.Info.MWTDB;
MWTDB = parseToleranceName(MWTDB);
D = MWTSet.Data_Plate;
D = innerjoin(MWTSet.Info.MWTDB,D,'Keys','mwtid');


%% TOLERANCE INDEX
%% open report
msrlist = {'speed','curve'};
fid1 = fopen(sprintf('%s/Effect Ind.txt',pSave),'w');
TAll = table;
posthocname = 'bonferroni';
pv = 0.05;
for msri = 1:numel(msrlist)
    %% add
    fprintf(fid1,'-- %s --\n',msrlist{msri});
    
    %% get 0mM-0mM control per exp mean
    a = MWTSet.Data_Exp.(msrlist{msri});
    i = ismember(a.predose,'0mM') & ismember(a.postdose,'0mM');
    a(~i,:) = [];
    % get control groups
    ag = strjoinrows([a.expname a.strain]);
    adata = a.mean;
    % create table
    RefT = table;
    RefT.group = ag;
    a = regexpcellout(RefT.group,' ','split');
    RefT.strain = a(:,2);
    RefT.a = adata;
    
    %% calculate tolerance
    % get 0mM 200mM exp mean
    c = MWTSet.Data_Exp.(msrlist{msri});
    i = ismember(c.predose,'0mM') & ismember(c.postdose,'200mM');
    c(~i,:) = [];
    % get control groups
    t = table; 
    t.group = strjoinrows([c.expname c.strain]);

    t.c = c.mean;
    RefT = innerjoin(RefT,t);
    % calculate A = c-a
    RefT.A = RefT.c - RefT.a;
    
    % get 200mM 200mM data
    d = D(ismember(D.postdose,'200mM') & ismember(D.predose,'200mM'),:);
    dg = strjoinrows([D.expname D.strain]);
    ddata = D.(msrlist{msri});
    t = table;
    t.group = dg;
    t.d = ddata;
    DT = innerjoin(t,RefT);
    % calculate tolerance
    DT.B = DT.d-DT.a;
    DT.TI = (DT.A-DT.B)./DT.A;
    
    %% calculate Withdrawal
    % determine k
    i = ismember(RefT.strain,'N2');
    m = RefT.c(i) > RefT.a(i);
    if sum(m) == sum(i)
        k = 1;
    elseif sum(~m) == sum(i)
        k = -1;
    else
        error('effect direction unstable in N2');
    end
    % get 200-0mM data
    d = D(ismember(D.postdose,'0mM') & ismember(D.predose,'200mM'),:);
    dg = strjoinrows([D.expname D.strain]);
    ddata = D.(msrlist{msri});
    t = table;
    t.group = dg;
    t.b = ddata;
    DW = innerjoin(t,RefT);
    % calculate withdrawal
    DW.WI = ((DW.b- DW.a)./DW.a).*k;
    
    
    %% export stats summary table
    % calculate stats summary
    T = grpstatsTable(DT.TI, DT.strain);
    T1 = table;
    T1.Index = repmat({'TI'},size(T,1),1);
    T1 = [T1 T];
    T = grpstatsTable(DW.WI, DW.strain);
    T2 = table;
    T2.Index = repmat({'WI'},size(T,1),1);
    T2 = [T2 T];
    T3 = [T1;T2];
    % add to all table
    T = table;
    T.msr = repmat(msrlist(msri),size(T3,1),1);
    T3 = [T T3];
    TAll = [TAll;T3];
    
    %% export anova
    [atext,phtext] = anova1_std(DT.TI,DT.strain);
    fprintf(fid1,'Tolerance Index:\nANOVA(strain): %s\n',atext);
    fprintf(fid1,'posthoc(%s) a=%.2f\n',posthocname,pv);
    fprintf(fid1,'%s\n\n',phtext);
    
    [atext,phtext] = anova1_std(DW.WI,DW.strain);
    fprintf(fid1,'Withdrawal Index:\nANOVA(strain): %s\n',atext);
    fprintf(fid1,'posthoc(%s) a=%.2f\n',posthocname,pv);
    fprintf(fid1,'%s\n\n',phtext);
    
end

fclose(fid1);
cd(pSave);
writetable(TAll,'Effect Index.csv');


return



















