function MWTSet = Stats_Output_ANOVA_GraphTable(Data,MWTSet,msrlist,component,pSave)
% create output

TAL = table; % table output
fid = fopen(fullfile(pSave,sprintf('ANOVA %s.txt',component)),'w'); % text output
for msri = 1:numel(msrlist)
    msr = msrlist{msri}; % get measure
    d = Data.(msr); % get data
    [txt,anovastats,multstats,T] = anova1_std_v2(d,Data.groupname,'plate');
    anovastats.pvalue = multstats.pValue;
    anovastats.descriptive = T;
    % add to central data
    MWTSet.Stats.(msr).(component).anova = anovastats;
    % text output
    fprintf(fid,'### %s ###\n%s\n\n\n',msr,txt);
    % table output
    n = size(T,1);
    msrt = cell2table(repmat({msr},n,1),'VariableNames',{'msr'});
    T = [msrt T];
    TAL  = [TAL;T];
end
fclose(fid); % close text output
writetable(TAL,fullfile(pSave,sprintf('Descriptive %s.csv',component)));
