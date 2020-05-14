function [StatOut,textfile] = RType_stats(T,varargin)

%% default
rmfactorName = 't';
factor1Name = 'gname';
pvlimit = 0.001;
alphanum = 0.05;
StatOut = [];
%%
vararginProcessor

%% prep for rmanova -------------------------
[TM,dvname,ivname, dmTable,~,dvvalues] = rmanova_extract_factor(T); % split gname into factors

%% start saving
textfile = '';

%% generate N
N = tabulateTable(TM.gname);
N = sortN2first(N.var, N);
a = strjoin(N.var,', ');
b = strjoin(num2cellstr(cell2mat(N.n))',', ');
textfile = sprintf('%sN(%s) = %s\n',textfile,a,b);

%%
% g = TM.gname;
% for coln = 1:size(dvvalues,2)
%    d = table2array(dvvalues(:,coln));
%    [n,m,sem] = grpstats(d,g,{'numel','mean','sem'});
%    
%    return
% end

% rm anova % --------------------
factorName = 'strain*dose';
astring = sprintf('%s-%s ~ %s',dvname{1},dvname{end},factorName);
rm = fitrm(TM,astring,'WithinDesign',dmTable);
StatOut.fitrm = rm;
t = ranova(rm); 
StatOut.ranova = t;
t = anovan_textresult(t,0, 'pvlimit',pvlimit);
textfile  = sprintf('%sRMANOVA(t:%s):\n%s\n',textfile,factorName,t);
% --------------------------------

% run pairwise % -------------------------------
astring = sprintf('%s-%s ~%s',dvname{1},dvname{end},factorName);
rm = fitrm(TM,astring,'WithinDesign',dmTable);
StatOut.fitrmg = rm;
t = multcompare(rm,factor1Name);
StatOut.mcomp_g = t;
t = multcompare_text2016b(t,'pvlimit',pvlimit,'alpha',alphanum);
textfile = sprintf('%s\nPosthoc(Tukey)curve by group:\n%s\n',textfile, t);
% -------------------------------------------

% comparison by taps -------------------------
% t = multcompare(rm,factor1Name,'By',rmfactorName);
% StatOut.mcomp_g_t = t;
% 
% t = multcompare_text2016b(t,'pvlimit',pvlimit,'alpha',alphanum);
% textfile = sprintf('%s\nPosthoc(Tukey)%s by %s:\n%s\n',textfile, factor1Name,rmfactorName,t);
textfile = sprintf('%s\nPosthoc(Tukey)%s by %s:\n',textfile, factor1Name,rmfactorName);
ssum = table;
for dvi =1:numel(dvname)
    dv = dvname{dvi};
    x= TM.(dv);
    group = TM.gname;
    [text,T,p,s,t,ST] = anova1_autoresults(x,group);
    [tt,gnames] = multcompare_convt22016v(s);
    
    result = multcompare_text2016b(tt,'grpnames',gnames,'prefix',[dv,'*']);
    textfile = sprintf('%s%s\n',textfile,result);
    
    t2 = table;
    a = str2num(regexprep(dv,'t',''));
    t2.t = repmat(a,size(tt,1),1);
    tt2 = [t2 tt];
    ssum = [ssum ; tt2];
end
StatOut.mcomp_g_t = ssum;
% --------------------------------------------

% comparison within group bewteen time --------
t = multcompare(rm,rmfactorName ,'By',factor1Name);
StatOut.mcomp_t_g = t;
t = multcompare_text2016b(t,'pvlimit',pvlimit,'alpha',alphanum);
textfile = sprintf('%s\nPosthoc(Tukey)%s by %s:\n%s\n',textfile, rmfactorName,factor1Name,t);
% --------------------------------------------


















