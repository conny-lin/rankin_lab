function [pstr, tv] = pairwisetime_comp(g1, g2, STG, rtime, alpha, pvlimit)


i = ismember(STG.gname_1,{g1,g2}) & ismember(STG.gname_2,{g1,g2});
a = STG(i,:);
t = a.t;
pv = a.pValue;
rval = pv(1:end) < alpha;
[pstr,~,tv] = pvaluestring2(rtime,rval,pv,alpha,pvlimit);