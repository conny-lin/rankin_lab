pHome = '/Users/connylin/Dropbox/RL/MWT_Analysis_ByExp/N2 10sISI/By Exp/N2 10sISI effect by exp';
pSave = fileparts(pHome);
asr = 'Initial';
[~,pE] = dircontent(pHome);
T = table;
for ei = 1:numel(pE)
    p = sprintf('%s/Dance_Glee_Showmance',pE{ei});
    [~,en] = fileparts(pE{ei});

    if isempty(dircontent(p,'Dance_Glee_Showmance.mat')) == 0;
        fn = sprintf('%s/Dance_Glee_Showmance/Dance_Glee_Showmance.mat',pE{ei});
        A = load(fn);
        
        t = table;
        t.expname = cellstr(en);
        a = A.MWTSet.Graph.(asr).RevFreq.Y;
        t.RevFreq0mM = a(1);
        if numel(a) == 2
            t.RevFreq400mM = a(2);
        else
            t.RevFreq400mM = nan;
        end
        t.RevFreqpvalue = A.MWTSet.Graph.(asr).RevFreq.ANOVA.Prob('Groups');
        t.RevFreqError = A.MWTSet.Graph.(asr).RevFreq.ANOVA.df('Error');
        a = A.MWTSet.Graph.(asr).RevSpeed.Y;
        t.RevSpeed0mM = a(1);
        if numel(a) == 2
            t.RevSpeed400mM = a(2);
        else
            t.RevSpeed400mM = nan;
        end
        t.RevSpeedpvalue = A.MWTSet.Graph.(asr).RevSpeed.ANOVA.Prob('Groups');
        t.RevSpeedError = A.MWTSet.Graph.(asr).RevSpeed.ANOVA.df('Error');
    else
        t.expname = cellstr(en);
        a = fieldnames(t);
        a(ismember(a,{'Properties','expname'})) = [];
        for x =1:numel(a)
            t.(a{x}) = nan;
        end
    end
    T = [T;t];
end
T.RevFreqD = T.RevFreq400mM - T.RevFreq0mM;
T.RevSpeedD = T.RevSpeed400mM - T.RevSpeed0mM;


writetable(T,sprintf('%s/RevFreq vs RevSpeed %s effects.csv',pSave,asr));

%% correlation
i = ~isnan(T.RevFreqD);
X = T.RevFreqD(i);
Y = T.RevSpeedD(i);
[rho,pvalue] = corr(X,Y);
df = numel(X)-2;
fid = fopen(sprintf('%s/RevFreq vs RevSpeed %s effects Correlation stat.txt',pSave,asr),'w');
if pvalue < 0.0001;
    fprintf('r(%d) = %.2f, p < 0.0001\n',df,rho);
    fprintf(fid,'r(%d) = %.2f, p < 0.0001\n',df,rho);
else
    fprintf('r(%d) = %.2f, p = %.3f\n',df,rho,pvalue);
    fprintf(fid,'r(%d) = %.2f, p = %.3f\n',df,rho,pvalue);

end
fclose(fid);
 
 
 

