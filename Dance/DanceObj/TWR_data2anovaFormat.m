function [DD,taps,A] = TWR_data2anovaFormat(MWTSet)

D = MWTSet.Raw;
mwtids = unique(D.mwtid);
taps = unique(D.tap);
A = nan(numel(mwtids),numel(taps));

for mi = 1:numel(mwtids)
    
   i = D.mwtid == mwtids(mi);
   d = D(i,:);
   d = sortrows(d,{'tap'});
   
   t = d.tap;
   y = d.RevFreq;
   
   [i,j] = ismember(taps,t);
   
   A(mi,j(i)) = y;
end

tap = A;
tap = array2table(tap);

% create legend
M = MWTSet.MWTDB(mwtids,:);
gn = M.groupname;
a = regexpcellout(gn,'_','split');
strain = a(:,1);
dose = a(:,2);
dose(cellfun(@isempty,dose)) = {'0mM'};

T = table;
T.groupname = gn;
T.strain = strain;
T.dose = dose;

DD = [T tap];