function T = parseRx(rx)


%% take out ones can be parsed

%% parse 
% NM1630_E3d12h400mM_R12h_T4d400mM
T = table;

%% test conditions
% rx = MWTDB.rx;
a = regexp(rx,'_','split');
[r,c] = cellfun(@size,a);
if size(a,1) > 1
    a = celltakeout(a);
else
    a = a{1};
end
% get last condition
b = cell(size(a,1),1);
for x = 1:size(a,1)
   b(x) = a(x,c(x));
end
T.age_test = cellfun(@str2num,regexpcellout(b,'\d(?=(d))','match'));
T.dose_test = cellfun(@str2num,regexpcellout(b,'\d{1,}(?=(mM))','match'));

%% recovery condition
T.rec_hr = nan(size(T,1),1);
a = regexpcellout(rx,'(?<=([_]R))\d{1,}(?=d)','match');
a(cellfun(@isempty,a)) = {'NaN'};
a = cellfun(@str2num,a);
T.rec_hr(~isnan(a)) = a(~isnan(a)).*24;
a = regexpcellout(rx,'(?<=([_]R))\d{1,}(?=h)','match');
a(cellfun(@isempty,a)) = {'NaN'};
a = cellfun(@str2num,a);
T.rec_hr(~isnan(a)) = a(~isnan(a));
a = regexpcellout(rx,'[_]Rvaries');
T.rec_hr(a==1) = Inf;

%% treatment condition (ignore second exposure)
% get rx condition
a = regexp(rx,'_','split');
[r,colN] = cellfun(@size,a);

a = regexpcellout(rx,'(?<=[_])R(\d{1,}\w{1,}|\w{1,})','match');
rxb4 = cell(size(rx));
for x = 1:numel(rx)
    rxb4(x) =  regexprep(rx(x),a(x),'');
end
rxb4 = regexprep(rxb4,'_','');

T.age_rx = cellfun(@str2num,regexpcellout(rxb4,'\d{1,}(?=d)','match'));
T.age_rx(colN==1) = NaN;

T.dose_rx = cellfun(@str2num,regexpcellout(rxb4,'\d{1,}(?=mM)','match'));
T.dose_rx(colN==1) = NaN;

%% duration rx
T.dur_rx = nan(size(T,1),1);
c = regexprep(rxb4,'mM','');
a = regexpcellout(c,'\d{1,}(?=h)','match');
a(cellfun(@isempty,a)) = {'NaN'};
a = cellfun(@str2num,a);
T.dur_rx(~isnan(a)) = a(~isnan(a)).*60;
a = regexpcellout(c,'\d{1,}(?=m)','match');
a(cellfun(@isempty,a)) = {'NaN'};
a = cellfun(@str2num,a);
T.dur_rx(~isnan(a)) = a(~isnan(a));
T.dur_rx(colN==1) = NaN;
