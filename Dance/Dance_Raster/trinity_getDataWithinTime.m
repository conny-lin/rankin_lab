function D = trinity_getDataWithinTime(DataAll,assayTime)

%% survey times
%  PROBLEM: bias can be -1, 1, or 0. 
% check if nan for bias
D = DataAll;
i = isnan(D.bias);
if sum(i) > 0
    D(i,:) = [];
    %     error('some bias record are NaN');
end
D(isnan(D.speed),:) = [];
% delete data outside of assay times
D(D.time < assayTime(1) | D.time >= assayTime(end),:) = []; 
