classdef rmANOVACurveReport
   properties
      datapath 
      curveData
      groupnames
      strain
      genotype
      n
      mean
      se
      posthoc
      manova
   end
   methods
        function curveData = get.curveData(obj)
            curveData = importfile_curveMANOVA(obj.datapath);
        end
        function n = get.n(obj)
            % find curve mean ++++++++
            a = obj.curveData{3};
            b = regexp(a,'(=)|[,]','split');
            c = regexprep(b(2:end),' ','');
            n = cellfun(@str2num,c);
            % ------------------------
        end
        function groupnames = get.groupnames(obj)
            % find group names ++++++++++++++
            a = obj.curveData{2};
            b = regexp(a,'(=)|[,]','split');
            c = regexprep(b(2:end),'^\s','');
            groupnames = c;
            % -------------------------------
        end
        function mean = get.mean(obj)
            % find curve mean ++++++++
            a = obj.curveData{4};
            b = regexp(a,'(=)|[,]','split');
            c = regexprep(b(2:end),' ','');
            mean = cellfun(@str2num,c);
            % ------------------------
        end
        function se = get.se(obj)
            % find curve mean ++++++++
            a = obj.curveData{5};
            b = regexp(a,'(=)|[,]','split');
            c = regexprep(b(2:end),' ','');
            se = cellfun(@str2num,c);
            % ------------------------
        end
        function manova = get.manova(obj)
            i = find(regexpcellout(obj.curveData,'Multifactorial')) +1;
            a = obj.curveData(i:i+2);
            % separate components
            b = regexpcellout(a,'(: )|(, )','split');
            T = table;
            T.factor = b(:,1);
            T.F = b(:,2);
            c = regexprep(b(:,3),' ','');
            T.ptxt = regexprep(c,'(<0.000)','<0.001');
            d = regexprep(c,'(p=)|(p<)','');
            d = regexprep(d,'(p=n.s.)','1');
            T.pvalue = cellfun(@str2num,d);
            % reconstruct sentence
            a = T.F;
            b =  T.factor;
            c = cellfunexpr(a,')');
            b = strjoinrows([c b],'');
            d = cellfun(@regexprep,a,c,b,'UniformOutput',0);
            d = strjoinrows([d T.ptxt],', ');
            T.txt = d;
            manova = T;
        end
        function posthoc = get.posthoc(obj)
            % find pairwise comparison +++++++
            i = find(regexpcellout(obj.curveData,'*** Posthoc tukey-kramer ***')) +1;
            a = obj.curveData(i:end);
            b = regexpcellout(a,'( x )|(,)','split');
            T = table;
            T.g1 = b(:,1);
            T.g2 = b(:,2);
            T.p = regexprep(b(:,3),'^\s','');
            T.pvalue = convertPvalueText2Num(T.p);
            posthoc = T;
            % -------------------------------
        end
        function strain = get.strain(obj)
            g = obj.groupnames;
            a = regexpcellout(g','^[A-Z]{1,}\d{1,}','match');
            b = unique(a);
            c(regexpcellout(b,'N2')) = [];
%             i = ~regexpcellout(g','N2');
            if isempty(c)
                strain = 'N2';
            elseif numel(c)==1
                strain = c;
            else
                error('more than one strain');
            end
%             strain = char(unique(regexprep(gs','( 0mM)|( 400mM)','')));
        end
        function pct_change = percentChange(obj)
            m = obj.mean;
            g = obj.groupnames;
            % calculate precent change: wildtype
            a = m(ismember(g,'N2 0mM'));
            b = m(ismember(g,'N2 400mM'));
            pct(1) = (((b-a)/a)*100);
            % calculate percent change: strain
            i = ~regexpcellout(g','N2');
            gs = g(i);
            ms = m(i);
            i = regexpcellout(gs','400mM');
            b = ms(i);
            a = ms(~i);
            pct(2) = (((b-a)/a)*100);
            % create gn
            gn = {'N2';obj.strain};
            % create table
            pct_change = array2table(pct','VariableNames',{'pct'},'RowNames',gn);
            % get pairwise comparison stats
            i = ismember(obj.posthoc.g1,'N2 0mM') & ismember(obj.posthoc.g2,'N2 400mM');
            if sum(i) ==1
            else
                i = ismember(obj.posthoc.g1,'N2 400mM') & ismember(obj.posthoc.g2,'N2 0mM');
            end
            p = obj.posthoc.p(i);
            i = ismember(obj.posthoc.g1,[obj.strain,' 0mM']) & ismember(obj.posthoc.g2,[obj.strain,' 400mM']);
            if sum(i) ==1
            else
                i = ismember(obj.posthoc.g1,[obj.strain,' 400mM']) & ismember(obj.posthoc.g2,[obj.strain,' 0mM']);
            end
            p(2) = obj.posthoc.p(i);
            pct_change.p = p';
            pct_change.pvalue = convertPvalueText2Num(pct_change.p);

        end
        function genotype = get.genotype(obj)
            s = obj.strain;
            strainNames = DanceM_load_strainInfo(s);
            i = ismember(strainNames.strain, s);
            genotype = strainNames.genotype{i};
        end
   end
    
    
end





























