function Tri = taptime_align_timed(Tri,timename)

%% find assay time alignment -----------------------------------
A = Tri(Tri.tap~=0,{'mwtid','time'});
A = unique(A,'rows');
% if multiple tap time found, use the first tap time
a = tabulate(A.mwtid);
i = a(:,2)>1;
if sum(i)>0
    mwtid = a(i,1);
    for mi = 1:numel(mwtid)
        t = min(A.time(A.mwtid==mwtid(mi))); % find first tap time
        A(A.mwtid==mwtid(mi) & A.time~=t,:) = []; % delete later times
    end
end
% record tap time
TapTime = A;

%% translate time to time away from taptime


B = Tri(:,{'mwtid','time'});
TapTime.taptime = TapTime.time;
TapTime.time = [];
A = innerjoin(B,TapTime);
A.td = A.time - A.taptime;
Tri.(timename) = A.td;
