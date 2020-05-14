classdef CurveStats_v1707
    properties
        MWTDB
        data
        mwtid
        groupnames
        expnames
        dscByPlate
        dscByGroup
        dscstat
        meanCtrl
        meanExp
        pctByPlate
        pctdscstat
        tstat
        statstr
    end
    methods
%         function mwtid = get.mwtid(obj)
%             mwtid = obj.Data.mwtid;
%         end
%         function data = get.data(obj)
%             data = obj.Data.data;
%         end
        function groupnames = get.groupnames(obj)
            a = obj.MWTDB(:,{'mwtid','groupname'});
            b = obj.Data(:,{'mwtid'});
            c = innerjoin(a,b);
            groupnames = c.groupname;
        end
        function expnames = get.expnames(obj)
            a = obj.MWTDB(:,{'mwtid','expname'});
            b = obj.Data(:,{'mwtid'});
            c = innerjoin(a,b);
            expnames = c.expname;
        end
        function dscByPlate = get.dscByPlate(obj)
            x = obj.data;
            g = obj.mwtid;
            M = obj.MWTDB;
            % make table
            S = statsBasicG(x,g,'mwtid');
            % get headers
            a = M(:,{'mwtid','expname','groupname','strain','rx'});
            b = S(:,{'mwtid'});
            L = innerjoin(b,a);
%             S.exp_gn = strjoinrows(table2cell(M(S.mwtid,{'expname','groupname'})),'*');
            dscByPlate = innerjoin(L,S);
        end
        function dscByGroup = get.dscByGroup(obj)
            A = obj.dscByPlate;
            g = A.groupname;
            x = A.mean;
            dscByGroup = statsBasicG(x,g,'groupname');
        end
        function meanCtrl = get.meanCtrl(obj)
            S = obj.dscByPlate;
            C = S(ismember(S.rx,'NA'),:);
            C.ctrl_name = strjoinrows([C.expname C.strain],'*');
            SC = statsBasicG(C.mean,C.ctrl_name,'ctrl_name');
            CL = unique(C(:,{'groupname','strain','rx','expname','ctrl_name'}));
            SC = innerjoin(CL,SC);
            SC.ctrl = SC.mean;
            SC.mean = [];
            meanCtrl = SC;
        end
        function meanExp = get.meanExp(obj)
            S = obj.dscByPlate;
            E = S(~ismember(S.rx,'NA'),:);
            E.ctrl_name = strjoinrows([E.expname E.strain],'*');
            meanExp = E;
        end
        function pctByPlate = get.pctByPlate(obj)
            E = obj.meanExp;
            C = obj.meanCtrl(:,{'ctrl_name','ctrl'});
            
            % divide against control
            A = innerjoin(E,C);
            A.pct = (A.mean./A.ctrl)-1;
            pctByPlate = A;
        end
        function [anovatext,T] = anova(obj)
            A = obj.pctByPlate;
            x = A.pct;
            g = A.groupname;
            if numel(x) ~= numel(g); error('check'); end
            if numel(unique(g))>1
                [anovatext,~,T] = anova1_std(x,g,'Tukey');
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
        end
        function tstat = get.tstat(obj)
            A = obj.pctByPlate;
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
            tstat = T;
        end
        function dscstat = get.dscstat(obj)
           [~,T] = anova(obj);
           dscstat = T;
        end
        function pctdscstat = get.pctdscstat(obj)
            A = obj.pctByPlate;
            x = A.pct;
            g = A.groupname;
            pctdscstat = statsBasicG(x,g,'groupname');
        end
        function statstr = get.statstr(obj)
            % descriptive
            A = obj.dscstat;
            g = char(strjoinrows(A.gnames',', '));
            n = char(strjoin(num2cellstr(A.N),', '));
            m = char(strjoin(num2cellstr(A.mean),', '));
            s = char(strjoin(num2cellstr(A.SE),', '));
            dsctxt = sprintf('G = %s\nN = %s\nmean = %s\nSE = %s',g,n,m,s);
            % anova
            [anovatext,~] = anova(obj);
            % ttest
            a = obj.tstat;
            b = strjoinrows([a.gn a.tstat a.pvs],', ');
            strttest = print_cellstring(b);
            % str
            statstr = sprintf('*** Descriptive ***\n%s\n\n*** ANOVA ***\n%s\n\n*** ttest 0 ***\n%s\n',dsctxt,anovatext,strttest);
        end
    end
    
end