function [DataGood,DataBad] = validate_trv_tap(Data)
%% CHECK TAP CONSISTENCY (r20151126)


%% get tap size
[r,c] = cellfun(@size,Data.data,'UniformOutput',0);
rn = cell2mat(r);

%% get exp name tap number
DB = parseMWTinfo(Data.mwtpath);
tapNExpected = DB.tapN;

%% find trv with inconsistent taps
ibad = rn~=tapNExpected;


%% remove bad taps
DataBad = Data(ibad,:);
DataGood = Data(~ibad,:);

%% report
pMWTBadTap = Data.mwtpath(ibad);
p = pMWTBadTap;
if isempty(p) == 0
d = parseMWTinfo(p);
fprintf('\nPlates with bad taps:\n');
tabulate(d.groupname)
end

p = DataGood.mwtpath;
d = parseMWTinfo(p);
fprintf('\nPlates with validated taps:\n');
tabulate(d.groupname)




