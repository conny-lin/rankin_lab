function [pstr,pvt,tv] = pvaluestring2(rtime,rval,pv,alpha,pvlimit)

%%
    a = [rtime rval pv];
    a((a(:,2) == 0),:) = [];
    a(a(:,3) > alpha,:) = [];
    tv = a(:,1);
    pvt = a(:,3);
    
    
    
    %%
    pstr = '';
    if isempty(pvt)
        pstr = 'p=n.s.';
    else
        if all(pvt < pvlimit)
           pstr = ['all, p<',num2str(pvlimit)];
        else
           n = size(pvt);
           str = arrayfun(@print_pvalue,pvt,...
               repmat(pvlimit,n),...
               repmat(alpha,n),'UniformOutput',0);
           t = num2cellstr(a(:,1));
           pstr = strjoinrows([t str],'s, ');
           pstr = char(strjoinrows(pstr',', '));
           
        end
    end
    
end