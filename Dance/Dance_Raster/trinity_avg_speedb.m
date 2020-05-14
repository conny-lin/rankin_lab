function [plateSumm,wormidU] = trinity_avg_speedb(D)


%% average speed * bias per frametimes(ftU)
% PROBLEM: sometimes speed has value as large as 1.32 but bias is 0. What
% does that mean?
% calculate speed dir
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

