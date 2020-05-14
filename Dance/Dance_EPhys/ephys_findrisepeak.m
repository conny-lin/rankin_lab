function [ymax,ymaxtime] = ephys_findrisepeak(d,t,upperbaseline)

% statsname = 'risepeak';
% for gi = 1:size(DataG,2)
%     i = (find(DataG(gi).time(1,:)==0)+1):timepts(gi);
%     d = nanmean(DataG(gi).speedbm(:,i));
ymax = max(d);
if ymax<=upperbaseline
    ymax = NaN; 
end
if ~isnan(ymax)
    ymaxtime = t(d==ymax);
else
    ymaxtime=NaN;
end
%     DataG(gi).(statsname) = ymax;
%     DataG(gi).([statsname,'_time']) = ymaxtime;
% end