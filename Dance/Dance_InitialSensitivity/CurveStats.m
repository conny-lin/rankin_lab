classdef CurveStats
    properties
        mwtid
        curve
        MWTDB
        groupnames
        statsByPlate
        statsByGroup
        pctEtohByPlate
        tstat
        statstr
    end
    methods
        function statsByPlate = get.statsByPlate(obj)
            x = obj.curve;
            g = obj.mwtid;
            M = obj.MWTDB;
            S = statsBasicG(x,g,'mwtid');
            S.exp_gn = strjoinrows(table2cell(M(S.mwtid,{'expname','groupname'})),'*');
            statsByPlate = S;
        end
        function pctEtohByPlate = get.pctEtohByPlate(obj)
            S = obj.statsByPlate;
            M = obj.MWTDB;
            A = obj.statsByPlate;

            x = S.mean;
            g = S.exp_gn;
            S1 = statsBasicG(x,g,'exp_gn');
            a = regexpcellout(S1.exp_gn,'*','split');
            S1.expname = a(:,1);
            S1.groupname = a(:,2);
            S1(regexpcellout(S1.groupname,'400mM'),:)= [];
            S0mM = S1;
            
            A = S;
            A(~regexpcellout(A.exp_gn,'400mM'),:) = [];
            A.exp_gn = regexprep(A.exp_gn,'_400mM','');

            [i,j] = ismember(A.exp_gn,S0mM.exp_gn);
            A.ctrlvalue = nan(size(A,1),1);
            A.ctrlvalue(i)= S0mM.mean(j(i));

            A.pct_ctrl = A.mean./A.ctrlvalue;
            B = table;
            B.groupname = M.groupname(A.mwtid);
            A(:,{'exp_gn','ctrlvalue'}) = [];
            pctEtohByPlate = [B A];
        end
        function [anovatext,T] = anova(obj)
            A = obj.pctEtohByPlate;
            x = A.pct_ctrl;
            g = A.groupname;
            if numel(x) ~= numel(g); error('check'); end
            if numel(unique(g))>1
                [anovatext,~,T] = anova1_std(x,g,'Tukey');
            else
                anovatext = 'no comparison made b/c only 1 group';
                % create desc stats
                T = statsBasic(A.pct_ctrl,'outputtype','struct');
                % correct naming convention
                T.gnames = unique(g);
                T.SE = T.se;
                T.N = T.n; 
                T.n = [];
            end
        end
        function tstat = get.tstat(obj)
            A = obj.pctEtohByPlate;
            gn = unique(A.groupname);
            gn = sortN2first(gn,gn);
            T = table;
            T.gn = gn;
            T.tstat = T.gn;
            T.pv = nan(size(T.gn,1),1);
            for gi =1:numel(gn)
                a = A.pct_ctrl(ismember(A.groupname,gn{gi}));
                [h,p,ci,stats] = ttest(a,0);
                T.tstat{gi} = sprintf('t(%d)=%.3f',stats.df,stats.tstat);   
                T.pv(gi) = p;
            end
            tstat = T;
        end
        function statstr = get.statstr(obj)
            % anova
            [anovatext,~] = anova(obj);
            % ttest
            a = obj.tstat;
            a.pv = cellfun(@print_pvalue,num2cell(a.pv),'UniformOutput',0);
            b = strjoinrows([a.gn a.tstat a.pv],', ');
            strttest = print_cellstring(b);
            % str
            statstr = sprintf('*** ANOVA ***\n%s\n\n*** ttest 0 ***\n%s\n',anovatext,strttest);
        end
    end
    
end