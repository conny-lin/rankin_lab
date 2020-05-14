function wrmval = extract_valid_wormtime(time,assaytimeSet)
%% sample input
% t1 = taptime-assaytimePretap;
% t2 = taptime;
% assaytimeSet = [t1;t2];
% time = Import.time;
%%
wrmval = false(size(time,1),size(assaytimeSet,2));
for ti = 1:size(assaytimeSet,2)
    t1 = assaytimeSet(1,ti);
    t2 = assaytimeSet(2,ti);
    wrmval(time(:,1) <= t1 & time(:,2) >= t2,ti) = true;
end
