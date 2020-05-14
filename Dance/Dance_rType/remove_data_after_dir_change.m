function [BSC,baseline_dir,baseline_dir_col] = remove_data_after_dir_change(BSBias,BS)


%% get baseline direction
[r,c] = ind2sub(size(BSBias),find(~isnan(BSBias)));
a = [r c];
[m,g] = grpstats(a(:,2),a(:,1),{'min','gname'});
g = cellfun(@str2num,g);
a = [g m];
i = sub2ind(size(BSBias),g,m);
if numel(m) ~= size(BSBias,1)
    error('all bias nan');
else
    baseline_dir = BSBias(i);
    baseline_dir_col = m;
end
V = BS.*BSBias;

%% get rid of baseine data in diff direction from the last baseline dir
BSC = BS;
BSC(baseline_dir == 1 & V<0) = NaN;
BSC(baseline_dir == -1 & V>0) = NaN;
BSC(baseline_dir == 0 & V~=0) = NaN;


%% find when dir switched
[r,c] = ind2sub(size(BSC),find(isnan(BSC)));
a = [r c];
%% get rid of ones before boundry
[i,j] = ismember(a(:,1),g);
a(i,3) = m(j(i));
a(:,4) = a(:,2)-a(:,3);
a(a(:,4) <=0,:) = [];
%% find boundry
a = sortrows(a(:,1:2),[1 2]);
boundry = [1;diff(a(:,1))];
a(~boundry,:) = [];



%% remove data after response dir switched 
for i = 1:size(a,1)
    r = a(i,1);
    c = a(i,2);
    BSC(r,c+1:end) = NaN;
end

%% check if all became nan
a = isnan(BSC);
a = sum(a,2);
if any(a==size(BSC,2))
    error('some data became all nan');
end

end