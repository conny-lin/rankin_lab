% calculate rapid tolerance index
% group name translation: specific for rapid tolerance
msrlist = {'speed','curve'};

MWTDB = MWTSet.Info.MWTDB;
MWTDB = parseToleranceName(MWTDB);
D = MWTSet.Data_Plate;
D = innerjoin(MWTSet.Info.MWTDB,D,'Keys','mwtid');


% k preset for withdrawal
k_withdrawal = {'speed',-1;'curve',1};

%% calculate EI
for msri = 1:numel(msrlist)
    % get WI, TI
    EM = MWTSet.Data_Exp.(msrlist{msri});
    %% get exp mean
    expstrain = strjoinrows([EM.expname EM.strain]);
    expstrainu = unique(expstrain);
    cond = strjoinrows([EM.predose EM.postdose]);
    condu = unique(cond);
    A = nan(numel(expstrainu),numel(condu));
    for ei = 1:numel(expstrainu)
        i = ismember(expstrain, expstrainu{ei});
        ED = EM(i,:);
        [i,j] = ismember(condu,cond(i));
        A(ei,i) = ED.mean(j(i));
    end
    RefT = table;
    RefT.group = expstrainu;
    a = regexpcellout(RefT.group,' ','split');
    RefT.a = A(:,1);
    RefT.c = A(:,2);
    RefT.b = A(:,3);
    RefT.d = A(:,4);
    RefT.A = RefT.c - RefT.a;
    
    %% get 200mM 200mM data (d, by plate)
    d = D(ismember(D.postdose,'200mM') & ismember(D.predose,'200mM'),:);
    t = table;
    t.group = strjoinrows([d.expname d.strain]);
    t.plate = d.mwtname;
    t.strain = d.strain;
    t.d = d.(msrlist{msri});
    DT = innerjoin(t,RefT(:,{'group','a','c','A'}));
    % calculate tolerance
    DT.B = DT.d-DT.a;
    DT.TI = (DT.A-DT.B)./DT.A;
    writetable(DT,sprintf('%s/%s TI raw.csv',pSave,msrlist{msri}));
    
    %% calculate Withdrawal
    % get 200-0mM data
    d = D(ismember(D.postdose,'0mM') & ismember(D.predose,'200mM'),:);
    t = table;
    t.group = strjoinrows([d.expname d.strain]);
    t.plate = d.mwtname;
    t.strain = d.strain;
    t.b = d.(msrlist{msri});
    DW = innerjoin(t,RefT(:,{'group','a','c'}));
    % calculate withdrawal
    % determine k
    k = k_withdrawal{ismember(k_withdrawal(:,1),msrlist{msri}),2};
    DW.WI = ((DW.b- DW.a)./DW.a).*k;
    
    %% write
    writetable(DT,sprintf('%s/%s TI raw.csv',pSave,msrlist{msri}));
    writetable(DW,sprintf('%s/%s WI raw.csv',pSave,msrlist{msri}));
    %% save to MWTSET
    MWTSet.EI.(msrlist{msri}).TI = DT;
    MWTSet.EI.(msrlist{msri}).WI = DW;
end


%% stats
% settings
fid1 = fopen(sprintf('%s/Effect Ind.txt',pSave),'w');
TAll = table;
TAllExp = table;
posthocname = 'bonferroni';
pv = 0.05;
% stats type
stype = {'plate','exp'};
% get data
Data = MWTSet.EI;
msrlist = fieldnames(Data);

% calculation
for msri = 1:numel(msrlist)
    %% output start
    fprintf(fid1,'-- %s --\n',msrlist{msri});  
    eiu = fieldnames(Data.(msrlist{msri}));
    for ei =1:numel(eiu)
        %% get data
        einame = eiu{ei};
        EI = Data.(msrlist{msri}).(einame);
        for stypei = 1:numel(stype)
            stypename = stype{stypei};
            fprintf(fid1,'%s, N(%s):\n',einame, stype{stypei});
            switch stypename
                case 'plate'
                    T = table;
                    T.strain = EI.strain;
                    T.mean = EI.(einame);
                case 'exp'
                    % calculate stats summary
                    T = grpstatsTable(EI.(einame), EI.group);
                    a = regexpcellout(T.gnameu,' ','split');
                    T.strain = a(:,2);
            end
            % anova
            if numel(unique(T.strain)) > 1
                [atext,phtext] = anova1_std(T.mean, T.strain);
                fprintf(fid1,'ANOVA(strain): %s\n',atext);
                fprintf(fid1,'posthoc(%s) a=%.2f\n',posthocname,pv);
                fprintf(fid1,'%s\n',phtext);
            end
            % t test
            fprintf(fid1,'%s: t test(right tail)\n',einame);
            gnu = unique(T.strain);
            gnu = [gnu(ismember(gnu,'N2')) gnu(~ismember(gnu,'N2'))];
            for gi = 1:numel(gnu)
                d = T.mean(ismember(T.strain,gnu(gi)));
                [h,p,~,STATS] = ttest(d,0,'tail','right');
                str = statTxtReport('t',STATS.df,STATS.tstat,p);
                fprintf(fid1,'%s, %s\n',gnu{gi},str);
            end 
            fprintf(fid1,'\n'); % enter another line

            % output table
            switch stypename
                case 'plate'
                    % output table
                    T = grpstatsTable(EI.(einame), EI.strain);
                case 'exp'   
                    % table
                    T = grpstatsTable(T.mean, T.strain);
            end
            
            % organize output table
            T.Properties.VariableNames(1) = {'strain'};
            i = ismember(T.strain,'N2');
            T = T([find(i) find(~i)],:);
            T1 = table;
            T1.msr = repmat(msrlist(msri),size(T,1),1);
            T1.Index = repmat({einame},size(T,1),1);
            T1 = [T1 T];            
            
            % add to master table
            switch stypename
                case 'plate'
                    TAll = [TAll;T1];
                case 'exp'
                    TAllExp = [TAllExp;T1];
            end
        end
    end
end

% write data
fclose(fid1);
cd(pSave);
writetable(TAll,'Effect Index Nplate.csv');
writetable(TAllExp,'Effect Index Nexp.csv');





%% export graphic data
TSum = struct;
TSum.plate = TAll;
TSum.exp = TAllExp;
stype = fieldnames(TSum);
for si = 1:numel(stype)
    Data = TSum.(stype{si});
    su = unique(Data.strain);
    su = [su(ismember(su,'N2'));su(~ismember(su,'N2'))];
    a = strjoinrows([su repmat({'se'},numel(su),1)],'_');
    var = [su;a];
    prefix = uniqueCellrows(Data(:,{'msr','Index'}));
    msru = unique(prefix.msr);
    indu = unique(prefix.Index);
    A = nan(size(prefix,1),numel(var));
    for msri = 1:numel(msru)
        for indi = 1:numel(indu)
            i = ismember(Data.msr,msru(msri)) & ismember(Data.Index,indu(indi));
            D = Data(i,:);
            [i,j] = ismember(D.strain,su);
            a = [D.mean(j(i))' D.se(j(i))'];
            i = ismember(prefix.msr, msru(msri)) & ismember(prefix.Index,indu(indi));
            A(i,:) = a;
        end
    end
    T = [prefix array2table(A,'VariableNames',var)];
    cd(pSave);
    writetable(T,sprintf('Effect Index N%s graph.csv',stype{si}));   
end

















