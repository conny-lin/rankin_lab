function T = combineStatsTable(M,SE,groupname,varargin)

%% default

%% treat input
if size(M,1) == numel(groupname)
   M = M';
   SE = SE';
end

%% transform into csv
a = groupname;
b = repmat({'Mean'},numel(a),1);
names1 = strjoinrows([b a],'_');
b = repmat({'SE'},numel(a),1);
names2 = strjoinrows([b a],'_');

T = array2table([M SE],'VariableNames',[names1 names2]);