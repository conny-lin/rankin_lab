function MWTSet = Stats_Output_ttest_pctCtrl(MWTSet,msrlist,Data,ctrl,component,pSave,nlimit)

%% get control by experiment
i = ismember(Data.groupname,ctrl);
DataC = Data(i,:);
% get plate N and exclude N<1
TN = tabulate(DataC.expname); 
n = cell2mat(TN(:,2));
expV = TN(n>nlimit,1);
Data(~ismember(Data.expname,expV),:) = [];
EXPList = cell2table(expV,'VariableNames',{'expname'});

%% split control and exp groups
i = ismember(Data.groupname,ctrl);
% get control data
DataC = Data(i,:); 
% get exp group data
DataE = Data(~i,:);
DC = EXPList;
for msri = 1:numel(msrlist)
    msr = msrlist{msri};
    D = statsBasicG(DataC.(msr),DataC.expname,'expname');
    D.Properties.VariableNames(ismember(D.Properties.VariableNames,'mean')) = {msr};
    DC = innerjoin(DC,D(:,{'expname',msr}));
end

%% calculate percent

% set up outputs
Data_PCT = DataE(:,{'mwtid','groupname','expname'});
DataE_expname =  DataE(:,{'expname'});
fid = fopen(fullfile(pSave,sprintf('ttest %s pct.txt',component)),'w');
T1 = table;
for msri = 1:numel(msrlist)
    msr = msrlist{msri}; % get msr name
    fprintf(fid,'%s\n',msr);

    % join control data with exp data
    A = innerjoin(DataE_expname,DC(:,{'expname',msr}));
    dpct = DataE.(msr)./A.(msr);
    Data_PCT.(msr) = dpct;
    MWTSet.Stats.(msr).(component).pct_raw = Data_PCT(:,{'mwtid','groupname','expname',msr});

    % stats by group by exp
    g = strjoinrows([Data_PCT.expname Data_PCT.groupname],'*');
    S = statsBasicG(Data_PCT.(msr), g,'v');
    a = regexpcellout(S.v,'*','split');
    S.v = a(:,2);
    S2 = statsBasicG(S.mean, S.v,'groupname');
    leg = cell2table(repmat({msr},size(S2,1),1),'VariableNames',{'msr'});
    T1 = [T1;[leg S2]]; 

    % stats per group
    glist = unique(S.v);
    for gi = 1:numel(glist)
        gn = glist{gi};
        SG = S(ismember(S.v,gn),:);

        stats = statsBasic(SG.mean,'outputtype','struct');
        [text,p,~,~] = ttest_auto(SG.mean,1);
        stats.group{gi} = gn;
        stats.pvalue(gi) = p;
        stats.ttest{gi} = text;
        stats.descriptive{gi} = stats;
        % text output
        fprintf(fid,'%s, %s\n',gn,text);
    end
    MWTSet.Stats.(msr).(component).pct_ttest = stats;

end
% table output
writetable(T1,fullfile(pSave,fullfile(sprintf('Descriptive PCT %s.csv',component))));
fclose(fid); % text output
