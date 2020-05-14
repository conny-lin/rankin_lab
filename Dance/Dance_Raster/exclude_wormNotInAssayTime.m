function D = exclude_wormNotInAssayTime(D,assayTime)
%% all worms must exist throughout the period
widU = unique(D.wormid);
t1 = assayTime(1);
t2 = assayTime(end-1);
n = nan(numel(widU),1);
nV = (numel(assayTime)-1);
for wi = 1:numel(widU)
    i = D.wormid == widU(wi);
    if ~isempty(i)
        t = D.frametime(i);
    end
    n = numel(unique(t));
    if n~=nV
        D(i,:) = []; % delete
    end
end
