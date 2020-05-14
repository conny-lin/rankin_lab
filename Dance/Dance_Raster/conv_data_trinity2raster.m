function [T,Array,Tlegend] = conv_data_trinity2raster(D,var,iv)

%  iv = {'assayperiod','wormid'};

%% TRANSFORM
R = struct;
switch var
    case 'velocity'
        D.velocity = D.speed./D.bias;
        D.velocity(D.bias==0)=0;
        [T,Array,Tlegend] = transposeTable(D,'velocity','frame',iv);
    otherwise
        [T,Array,Tlegend] = transposeTable(D,var,'frame',iv);

end
return


%% translate frametime(ft) into assay times
ft = D.time;
for ti = 1:numel(assayStartTime)-1
%     if ti == numel(assayTime)
%         i = ft >= assayTime(ti);
%     else
        i = ft >= assayStartTime(ti) & ft < assayStartTime(ti+1);
%     end
    if sum(i) == 0
       [~,fn] = fileparts(pMWTp);
       plateSumm = [];
       warning('time point (%.2fs) in plate [%s] has no data, skip',assayStartTime(ti),fn); 
       return
    end
    ft(i) = assayStartTime(ti);
end
% validate:
% 1. all frame time must equal to one of the assay time
% 2. number of unique frame time must equal to (assay time -1)
if sum(ismember(ft,assayStartTime)) ~= numel(ft) ||...
        numel(unique(ft)) ~= numel(assayStartTime)-1
    error('frame time matching failed');
else
    D.frametime = ft;
end


%% all worms must exist throughout the period
widU = unique(D.wormid);
t1 = assayStartTime(1);
t2 = assayStartTime(end-1);
n = nan(numel(widU),1);
nV = (numel(assayStartTime)-1);
for wi = 1:numel(widU)
    i = D.wormid == widU(wi);
    t = D.frametime(i);
    n = numel(unique(t));
    if n~=nV
        D(i,:) = []; % delete
    end
end


%% average speed * bias per frametimes(ftU)
% PROBLEM: sometimes speed has value as large as 1.32 but bias is 0. What
% does that mean?
% calculate speed dir
% fprintf('-converting speed to velocity...');
D.speedDir = D.speed.*D.bias;
wormidU = unique(D.wormid);
Nwormid = numel(wormidU);
ftU = unique(D.frametime);
A = nan(Nwormid,numel(ftU));
for x = 1:numel(ftU)
    t = ftU(x);

    for wi = 1:numel(wormidU)
        d = D(D.frametime == t & D.wormid == wormidU(wi),:);
        if sum(diff(d.bias)) == 0 
        % if move within this frame all in one dir
        % calculate mean of all data
            A(wi,x) = mean(d.speedDir);
        elseif sum(diff(d.bias)) ~= 0 
        % if movement within this frame contains dir shift
            if sum(d.bias < 0) > 0
            % if reversal exists, only calculate reversal mean
                A(wi,x) = mean(d.speedDir(d.bias < 0));
            elseif sum(d.bias > 0) > 0
            % if only forward exists, only calculate forward mean
                A(wi,x) = mean(d.speedDir(d.bias > 0));
            elseif sum(d.bias == 0) > 0
            % if only bias == 0 (presumably pause), only calculate pase
            % speedDir, which is set to zero
                A(wi,x) = mean(d.speedDir(d.bias == 0));
            end
        end
    end
end
plateSumm = A;
% if any nan results, flag
if sum(any(isnan(A))) > 0
   error('some nan results, code to fix'); 
end
% 


%% output
% create output structure array
% A = struct;
% A.wormID = wormID;
% A.frametime = frameTimeRecord;
% A.data = plateSumm;
% A.pMWT = pMWTp;

% text output
if saveopt == 1
    dlmwrite(sprintf('%s/%ds_%ds_fint_%.1f_N%d_wormID.trinitySummary',pMWTp,start,finish,frameInt,size(wormID,1)),wormID);
    dlmwrite(sprintf('%s/%ds_%ds_fint_%.1f_N%d_frametime.trinitySummary',pMWTp,start,finish,frameInt,size(wormID,1)),frameTimeRecord)
    dlmwrite(sprintf('%s/%ds_%ds_fint_%.1f_N%d_rasterdata.trinitySummary',pMWTp,start,finish,frameInt,size(wormID,1)),plateSumm)
end

