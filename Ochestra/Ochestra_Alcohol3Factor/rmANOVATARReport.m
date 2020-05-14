classdef rmANOVATARReport
   properties
      datapath 
      txtcontent
      strain
      genotype
      anova
      posthoc_groups
      posthoc_tapbygroup
      posthoc_tapbygroupStr
   end
   methods
        function txtcontent = get.txtcontent(obj)
            txtcontent = importfile_curveMANOVA(obj.datapath);
        end
        function anova = get.anova(obj)
            D = obj.txtcontent;
            i = find(regexpcellout(D,'RMANOVA')) +2;
            j = find(regexpcellout(D,'Posthoc[(]Tukey[)]')) -1;
            fulltxt = D(i:j);
            % separate components
            b = regexpcellout(fulltxt,'(: )|(, )','split');
            ptxt = b(:,2);
            a = b(:,1);
            b = regexprep(a,'*tap','');
            c = regexpcellout(b,'([)])|( = )','split');
            rownames = c(:,2);
            T = table;
            T.factor = rownames;
            T.fulltxt = fulltxt;
            T.F = b;
            T.ptxt = ptxt;
            T.Properties.RowNames = rownames;
            d = regexprep(T.ptxt,'(p=)|(p<)','');
            d = regexprep(d,'(n.s.)','1');
            T.pvalue = cellfun(@str2num,d);
            anova = T;
        end
        function posthoc_groups = get.posthoc_groups(obj)
            D = obj.txtcontent;
            i = find(regexpcellout(D,'Posthoc[(]Tukey[)]'));
            j = i(1)+1;
            k = i(2)-1;
            fulltxt = D(j:k);

            b = regexpcellout(fulltxt,'(*)|(, )','split');
            T = table;
            T.g1 = b(:,1);
            T.g2 = b(:,2);
            T.ptxt = b(:,3);

            b = regexpcellout(fulltxt,'(, )','split');
            T.Properties.RowNames = b(:,1);
            T.pvalue = convertPvalueText2Num(T.ptxt);
            posthoc_groups = T;
        end
        function strain = get.strain(obj)
            g = obj.posthoc_groups.g1;
            g(regexpcellout(g,'N2')) = [];
            g = regexprep(g,'_400mM','');
            strain =unique(g);
        end
        function posthoc_tapbygroup = get.posthoc_tapbygroup(obj)
            D = obj.txtcontent;
            
            i = find(regexpcellout(D,'Posthoc[(]Tukey[)]'));
            fulltxt = D(i(2)+1:end);
            b = regexpcellout(fulltxt,'(*)|(, )','split');
            T = table;
            T.tap = str2numcell(regexprep(b(:,1),'tap',''));
            T.g1 = b(:,2);
            T.g2 = b(:,3);
            T.ptxt = b(:,4);
            r = strjoinrows(b(:,1:3),'*');
            T.Properties.RowNames = r;
            T.pvalue = convertPvalueText2Num(T.ptxt);            
            posthoc_tapbygroup = T;
        end
        function str = genPosthocTapbyGStr(obj,g1,g2,pvlim,pvsig)
            S = obj.posthoc_tapbygroup;
            i = ismember(S.g1,g1) & ismember(S.g2,g2);
            if sum(i) == 0
                i = ismember(S.g1,g2) & ismember(S.g2,g1);
            end
            S1 = S(i,:);
            pstr = pvaluestring(S1.tap,S1.pvalue,pvsig,pvlim);
            s = strjoin([{g1 g2}],' vs ');
            s = regexprep(s,'N2','wild-type');
            s = regexprep(s,'_',' ');
            s = regexprep(s,obj.strain,obj.genotype);
            str = [s ,', ',pstr];
            
        end
        function genotype = get.genotype(obj)
            strainNames = DanceM_load_strainInfo(obj.strain);
            genotype = strainNames.genotype;
        end

   end
    
    
end





























