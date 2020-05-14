        function pv = convertPvalueText2Num(ptxt)
            d = regexprep(ptxt,'(p=)|(p<)','');
            d = regexprep(d,'(n.s.)','1');
            pv = cellfun(@str2num,d);
        end