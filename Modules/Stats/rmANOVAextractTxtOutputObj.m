classdef rmANOVAextractTxtOutputObj
   properties
      datapath 
      txtcontent
      strain
      genotype
      groups
      anova
      posthoc_groups
      txt_t_identifier
      posthoc_tXgroup
%       posthoc_tXgroupStr
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
        function groups = get.groups(obj)
                groups = unique([obj.posthoc_groups.g1; obj.posthoc_groups.g2]);
 
        end
        function strain = get.strain(obj)
            g = obj.posthoc_groups.g1;
            g(regexpcellout(g,'N2')) = [];
            g = regexprep(g,'_400mM','');
            strain =unique(g);
        end
        function T = get.posthoc_tXgroup(obj)
            D = obj.txtcontent;
            txtsep = obj.txt_t_identifier;
            i = find(regexpcellout(D,txtsep));
            fulltxt = D(i(1)+1:end);
            b = regexpcellout(fulltxt,'(*)|(, )','split');
            T = table;
            T.t = str2numcell(regexpcellout(b(:,1),'\d{1,}','match'));
            T.g1 = b(:,2);
            T.g2 = b(:,3);
            T.g1xg2 = strjoinrows([T.g1,T.g2],'*');
            T.ptxt = b(:,4);
            r = strjoinrows(b(:,1:3),'*');
            T.Properties.RowNames = r;
            T.pvalue = convertPvalueText2Num(T.ptxt);
            T = sortrows(T,{'g1xg2','t'});        
        end
        function str = posthoc_tXgroup_pairStr(obj,g1,g2,pvlim,pvsig)
            S = obj.posthoc_tXgroup;
            i = ismember(S.g1,g1) & ismember(S.g2,g2);
            if sum(i) == 0
                i = ismember(S.g1,g2) & ismember(S.g2,g1);
            end
            S1 = S(i,:);
            pstr = pvaluestring(S1.t,S1.pvalue,pvsig,pvlim);
            s = strjoin([{g1 g2}],' vs ');
            s = regexprep(s,'N2','wild-type');
            s = regexprep(s,'_',' ');
            s = regexprep(s,obj.strain,obj.genotype);
            str = [s ,', ',pstr];
        end
        function str = posthoc_tXgroupWriteUp(obj,pvlim,pvsig)
            pairwise_group_pairs = pairwisecomp_getpairs(obj.groups);
            str = sprintf('');
            for x = 1:size(pairwise_group_pairs,1)
                g1 = pairwise_group_pairs{x,1};
                g2 = pairwise_group_pairs{x,2};
                t = posthoc_tXgroup_pairStr(obj,g1,g2,pvlim,pvsig);
                str = sprintf('%s\n%s',str,t);
            end
        end
        function genotype = get.genotype(obj)
            strainNames = DanceM_load_strainInfo(obj.strain);
            genotype = strainNames.genotype;
        end

   end
    
    
end





























