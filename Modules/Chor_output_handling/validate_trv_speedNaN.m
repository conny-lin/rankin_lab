function [DataGood,DataBad] = validate_trv_speedNaN(Data)


%% calculate speed
v = false(size(Data,1),1);
for x = 45:size(Data,1)
    d = Data.data{x};
    n = d.N_ForwardOrPause + d.N_Rev;
    a = d.N_Rev./n;
    % calculate speed = RevDis/RevDur
    a = d.RevDis./d.RevDur;
    % if all tap has value, record as valid
    v(x) = ~any(isnan(a));
    return
end

%% remove 
DataBad = Data(~v,:);
DataGood = Data(v,:);

