function assayTime = trinity_calAssaytime(D,start,finish,varargin)

%%
version = 'cleaned';
frameInt = 0.2; %(ed20161004)
vararginProcessor

switch version
    case 'uncleaned'
        %% calculate assay time
        % tap time is indicated as colume 7 == 1, time is at column 1
        % find tap time
        a = 0;
        t = [];
        for wi = 1:size(D,1)
            i = D{wi,2}(:,7) ==1;
            if  sum(i) > 0
                t = [t; D{wi,2}(i,1)];
            end
        end
        tapTime = unique(t);
        

    case 'cleaned'
        %%
        legend = {'wormid','time','tap','speed_bias'};
        if ~istable(D)
            D = array2table(D,'VariableNames',legend);
        end
        tapTime = unique(D.time(D.tap == true));

        
end

%%
if isempty(tapTime) == 0 % if tap is found
    % find taps within start and finish time
    t = tapTime;
    i = t > start & t < finish;
    t = t(i);
    if isempty(t) == 0 % if one or more tap is found, 
        % align to the first tap
        tstart = t(1);
        assayTime = [flip(tstart-frameInt:-frameInt:start)  tstart:frameInt:finish+frameInt];
    elseif sum(i) == 0 % if no tap found in the time frame, 
        % align to first tap
        df = tapTime(1) - floor(tapTime(1));
        tstart = start + df;
        assayTime = [flip(tstart-frameInt:-frameInt:start)  tstart:frameInt:finish+frameInt];
    end
else % if no tap is found
    assayTime = start:frameInt:finish;
end








