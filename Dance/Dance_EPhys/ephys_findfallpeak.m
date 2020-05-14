function [ymin,ymintime] = ephys_findfallpeak(d,t,lowerbaseline)


ymin = min(d);
if ymin>=lowerbaseline
    ymin = NaN; 
end
if ~isnan(ymin)
    ymintime = t(d==ymin);
else
    ymintime=NaN;
end
