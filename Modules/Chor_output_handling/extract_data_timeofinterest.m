function Data = extract_data_timeofinterest(time,Data,t)

% keep only time of interest
Data(time < t(1) | time > t(2),:) = [];