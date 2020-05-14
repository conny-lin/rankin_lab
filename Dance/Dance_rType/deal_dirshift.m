function S = deal_dirshift(S,timename,varargin)

vararginProcessor;

%% deal with dir shifts within frame
Raw = S.bias.timetable_raw;

idirshift = find(~ismember(S.bias.array,[1 0 -1]) & ~isnan(S.bias.array));
arraySize = size(S.bias.array);

for di = 1:numel(idirshift)
    iissue = idirshift(di);
    [r,c] = ind2sub(arraySize,iissue);
    t = S.bias.timetable.(timename)(r);
    wormid = char(S.bias.timetable.Properties.VariableNames(c));
    speed = table2array(S.speed.timetable(t,wormid));
    bias = table2array(S.bias.timetable(t,wormid));
    if ~isnan(speed) || ~isnan(bias)
        r1 = r-1;
        r2 = r+1;
        if r1 <1; r1 = 1; end
        if r2 > arraySize(1), r2 = arraySize(1); end

        i = ismember(Raw.ids, wormid) & ...
            Raw.(timename) >= S.bias.timetable.(timename)(r1) &...
            Raw.(timename) <= S.bias.timetable.(timename)(r2);
        R = Raw(i,:);
        if ~isempty(R)

            if any(R.bias == -1) % if any bias = reversal, use reversal mean
                k = i(R.bias == -1);
                s = nanmean(S.speed.timetable_raw.speed(k));
                S.bias.array(iissue) = -1;
            elseif ~any(R.bias == -1) && any(R.bias == 1) % if no reversal and has forward, use forward mean
                k = i(R.bias == 1);
                s = nanmean(S.speed.timetable_raw.speed(k));
                S.bias.array(iissue) = 1;
            elseif ~any(R.bias == -1) && ~any(R.bias == 1) && any(R.bias ==0) % if no reversal and has forward, use forward mean
                k = i(R.bias == 0);
                s = nanmean(S.speed.timetable_raw.speed(k));
                S.bias.array(iissue) = 0;
            else
                [r c]
                error('not accomodating')
            end
        else
            error('not accomodating')

        end
        S.speed.array(iissue) = s;
    end
    
   
end
%% check if successful
i = ~ismember(S.bias.array,[1 0 -1]) & ~isnan(S.bias.array);
if sum(sum(i)) >0
    error('still some ambiguous bias');
end

%%

S.speed.timetable = array2timetable(S.speed.array,'RowTimes',S.speed.timetable.(timename),'VariableNames',S.speed.timetable.Properties.VariableNames);
S.bias.timetable = array2timetable(S.bias.array,'RowTimes',S.bias.timetable.(timename),'VariableNames',S.bias.timetable.Properties.VariableNames);




