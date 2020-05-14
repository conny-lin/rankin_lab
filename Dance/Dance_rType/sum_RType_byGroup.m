function [accAvg, accSem, accN, tt, MWD] = sum_RType_byGroup(MWTDB, mwtidlist,A,NV)

%% check input
gnlist = unique(MWTDB.groupname);

[i,j] = ismember(mwtidlist, MWTDB.mwtid);
if sum(i) ~= numel(mwtidlist)
    error('some mwt has no information'); 
end
plateidu = j(i);
MWD = MWTDB(plateidu,:);


%%

accAvg = nan(numel(gnlist),size(NV,2));
accSem = accAvg;
accN = accAvg;
for gi = 1:numel(gnlist) % cycle through groups
    gn = gnlist{gi}; % get gname

    i = ismember(MWD.groupname,gn);

    if ~isempty(i)
%         a = A(i,:);
%         n = NV(i,:);
%         j = ~any(n < n_lowest,2);

        aa = A(i,:); % get data matching group name
        ap = nanmean(aa);        
        accAvg(gi,:) = ap;
        
        nn = size(aa,1);
        nn = sum(~isnan(aa));
        
        accN(gi,:) = nn;
        
        apse = nanstd(aa);
        apsem = apse./sqrt(nn-1);
        accSem(gi,:) = apsem;
    end
end

t = accAvg;
t = array2table(t);
tt = table;
tt.gname = gnlist;
tt = [tt t];