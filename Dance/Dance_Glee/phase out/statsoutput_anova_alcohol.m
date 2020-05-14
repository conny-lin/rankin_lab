function A = statsoutput_anova_alcohol(prefix,A,B,G,mwtnames,pSaveA,msr)
    %% create output
    S = table;
    S.mwtname = mwtnames;
    S.groupname = G;
    % break G into strain vs alcohol
    a = regexpcellout(G,'(\<\w+\d+(?=[_]\w*))|(\<\w+\d+\>)','match');
    b = regexpcellout(G,'(?<=\<\w{1,}\d{1,}_)\w*','match');
    b(cellfun(@isempty,b)) = {'0mM'};
    G2 = {a(:,1),b};  
    st = a(:,1);
    S.strain = st;
    S.dose = b;
    S.data = B;
    A.(msr).Raw = S;
    if numel(unique(S.groupname)) > 2
        % anovan
        [~,t,~] = anovan(B,G2,'varnames',{'strain','dose'},'model','full','nested',[0 0;1 0],'display','off');
        [r,T] = anovan_textresult(t);
        writetable(T,sprintf('%s/%s %s ANOVA.csv',pSaveA,msr,prefix),'Delimiter',',');
        A.(msr).ANOVAN = T;
        % export anovan text output
        fid = fopen(sprintf('%s/%s %s ANOVA.txt',pSaveA,msr,prefix),'w');
        fprintf(fid,'Multifactorial nested ANOVA\n');
        for x = 1:numel(r); fprintf(fid,'%s\n',r{x}); end
        fclose(fid);
        % posthoc
        [~,~,stats] = anova1(B,G,'off');
        [c,~,~,gnames] = multcompare(stats,'ctype','bonferroni','display','off');
        [T,r] = multcompare_pairinterpretation(c,gnames,'nsshow',0);
        A.(msr).posthoc = T;
        A.(msr).posthoctest = 'bonferroni';
        writetable(T,sprintf('%s/%s %s posthoc bonferroni.csv',pSaveA,msr,prefix),'Delimiter',',');
        % export posthoc text output
        fid = fopen(sprintf('%s/%s %s posthoc bonferroni.txt',pSaveA,msr,prefix),'w');
        for x = 1:numel(r); fprintf(fid,'%s\n',r{x}); end
        fclose(fid);
    else 
        % anova
        [~,t,~] = anova1(B,G,'off');
        [r,T] = anovan_textresult(t);
        writetable(T,sprintf('%s/%s %s ANOVA.csv',pSaveA,msr,prefix),'Delimiter',',');
        A.(msr).ANOVA = T;
        % export anovan text output
        fid = fopen(sprintf('%s/%s %s ANOVA.txt',pSaveA,msr,prefix),'w');
        fprintf(fid,'ANOVA %s\n',msr);
        for x = 1:numel(r); fprintf(fid,'%s\n',r{x}); end
        fclose(fid);
    end
end
