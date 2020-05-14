function [DataGood,DataBad] = validate_trv_freqNaN(Data)


%%
v = false(size(Data,1),1);
for x = 1:size(Data,1)
    d = Data.data{x};
    n = d.N_ForwardOrPause + d.N_Rev;
    a = d.N_Rev./n;
    i = isnan(a);
    v(x) = sum(i) == 0;
end

%% remove 
DataBad = Data(~v,:);
DataGood = Data(v,:);

