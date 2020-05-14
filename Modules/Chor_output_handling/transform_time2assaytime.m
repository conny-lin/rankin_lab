function timea = transform_time2assaytime(time,timeAssay)

timea = nan(size(time));
for ti = 1:numel(timeAssay)-1
    t1 = timeAssay(ti);
    t2 = timeAssay(ti+1);
    i = time > t1 & time <= t2;
    timea(i) = t2;
end