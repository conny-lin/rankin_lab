function [pv,tr] = get_pvalue(ST,tmaxi,gname)

i = ismember(ST.groupname,gname) & ST.t_1 == 0 & ismember(ST.t_2,1:tmaxi);
pv = ST.pValue(i);
tr = ST.t_2(i);
    