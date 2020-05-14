function str = taptime_pvalue(tv,pstr)

%%
if numel(tv)>1
    tvs = tv(1):0.1:tv(end);
    n = numel(tvs);
    if numel(tv) < n
       a = num2cellstr(tv')';
       b = repmat({'s'},numel(a),1);
       a = strjoinrows([a b],'');
       str = strjoin(a',',');
       str = sprintf('%s after the tap (%s)',str,pstr);
   else
       t1 = tv(1);
       t2 = tv(end);
       str = sprintf('%.1f-%.1fs after the tap (%s)',t1,t2,pstr);

    end
elseif numel(tv)==1
    str = sprintf('%.1fs after the tap (%s)',tv,pstr);
elseif isempty(tv)
    str ='';
end