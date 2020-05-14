function [timestd,tapN] = syncdata2tap(time,tap,taptimestd)
if nargin<=2
    taptimestd = 0;
end
timestd = nan(size(time));

%%
taptime = time(tap==1);
tapN = numel(taptime);
if tapN>1;
    warning('code not able to accomodate multiple taps');
elseif tapN==0
    warning('no tap found');
else
    timestd = time-taptime+taptimestd;
end

