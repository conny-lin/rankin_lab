function StatOut = rmanova_TWR(T,varargin)

%%
rmfname = 'tap';
factor1Name = 'groupname';
a = T.Properties.VariableNames';
a = regexpcellout(a,['(?<=',rmfname,')\d{1,}'],'match');
a(cellfun(@isempty,a)) = [];
a = cell2mat(cellfun(@str2num,a,'UniformOutput',0));

rmtable = array2table(a,'VariableNames',{rmfname});

t1 = rmtable.(rmfname)(1);
t2 = rmtable.(rmfname)(end);


pvlimit = 0.001;
alphanum = 0.05;

vararginProcessor;

StatOut = struct;

%%


%% anova
factorName = 'strain*dose';
astring = sprintf('%s%d-%s%d~%s',rmfname,t1,rmfname,t2,factorName);
rm = fitrm(T,astring,'WithinDesign',rmtable);
t = ranova(rm); 
StatOut.ranova = t;
% t = anovan_textresult(t,0, 'pvlimit',pvlimit);
% fprintf(fid1,'RMANOVA(%s,t:%s):\n%s\n',timeName,factorName,t);

% run pairwise
astring = sprintf('%s%d-%s%d~%s',rmfname,t1,rmfname,t2,factor1Name);
rm = fitrm(T,astring,'WithinDesign',rmtable);
t = multcompare(rm,factor1Name);
StatOut.mcomp_g = t;
% t = multcompare_text2016b(t,'pvlimit',pvlimit,'alpha',alphanum);
% fprintf(fid1,'\nPosthoc(Tukey)curve by group:\n%s\n',t);


% comparison by taps
t = multcompare(rm,factor1Name,'By',rmfname);
StatOut.mcomp_g_t = t;
t = multcompare_text2016b(t,'pvlimit',pvlimit,'alpha',alphanum);
% fprintf(fid1,'\nPosthoc(Tukey)%s by %s:\n%s\n',factor1Name,rmfactorName,t);

% comparison within group bewteen time
t = multcompare(rm,rmfname ,'By',factor1Name);
StatOut.mcomp_t_g = t;
% t = multcompare_text2016b(t,'pvlimit',pvlimit,'alpha',alphanum);
% fprintf(fid1,'\nPosthoc(Tukey)%s by %s:\n%s\n',rmfactorName,factor1Name,t);

% fprintf(fid1,'\n\n');
% ----------------------





