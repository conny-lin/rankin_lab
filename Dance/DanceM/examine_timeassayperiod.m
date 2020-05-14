function [A,tval,a,b,period1,period2]= examine_timeassayperiod(t1,t2,ti,tf)

% ti = floor(TimeSet.atstart); % start time
% tf = ceil(TimeSet.atend); % end time
% 
% t1 = TrinityInfo.t1;
% t2 = TrinityInfo.t2;


tim = repmat(ti,size(t1,1),1);
period1 = (numel(ti)+1) - sum(t1 <= tim,2); % this get you the first period of assay

tfm = repmat(tf,size(t2,1),1);
period2 = sum(t2 >= tfm,2); % this get you the last period of assay

tval = true(size(t1,1),1);
% tval(period1 <= 0 | (period2 ==0 & period1~=1)) = false;
tval(period1 <= 0 | period2 ==0) = false;

period1(period1>2) = period1(period1>2)-1;

a = [(1:numel(t1))' period1 period2 t1 t2 tval];
b = a(~tval,:);

A = array2table(a,'VariableNames',{'id','period1','period2','t1','t2','tval'});

a(~tval,:) = [];
