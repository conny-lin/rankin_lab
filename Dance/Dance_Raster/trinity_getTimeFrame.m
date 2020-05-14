function D = trinity_getTimeFrame(D,assayTime,pMWTp)


%% translate frametime(ft) into assay times
ft = D.time;
for ti = 1:numel(assayTime)-1

    i = ft >= assayTime(ti) & ft < assayTime(ti+1);
    if sum(i) == 0

       [~,fn] = fileparts(pMWTp);
       warning('time point (%.2fs) in plate [%s] has no data, skip',assayTime(ti),fn); 
       D = table;

    else
       
        ft(i) = assayTime(ti);
        % validate:
        % 1. all frame time must equal to one of the assay time
        % 2. number of unique frame time must equal to (assay time -1)
        if sum(ismember(ft,assayTime)) ~= numel(ft) ||...
                numel(unique(ft)) ~= numel(assayTime)-1
            warning('frame time matching failed');
        else
            D.frametime = ft;

        end

    end
end   
