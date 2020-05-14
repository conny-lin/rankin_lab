function [TrinityInfo,TimeSet] = get_trinityIndInfo(TrinityInfo,TimeSet)



% get trinity files interested within time
ti = floor(TimeSet.atstart); % start time
tf = ceil(TimeSet.atend); % end time
t1 = TrinityInfo.t1;
t2 = TrinityInfo.t2;
A = examine_timeassayperiod(t1,t2,ti,tf);
TrinityInfo = [TrinityInfo A(:,{'period1','period2','tval'})];
TrinityInfo(~A.tval,:) = []; % trim files
