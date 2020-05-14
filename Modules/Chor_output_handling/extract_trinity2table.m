function D = extract_trinity2table(TrinityData,var)


%% get legend
Legend = chormaster5('TrinityOnly','','legendonly',true);
Legend = Legend.trinity;
Legend = [{'wormid'} ;Legend];
if nargin<2
    var = Legend;
end

%% calculate assay time
% tap time is indicated as colume 7 == 1, time is at column 1

[i,j] = ismember(var,Legend);
if sum(i)~=numel(i)
    error('no required var');
end
ivar = j(i);

%% generate worm id
D = TrinityData;
wormN = cell2mat(cellfun(@size,D(:,2),'UniformOutput',0)); % get size
wormNc = cell(size(D,1),1);
for i= 1:size(D,1)
    a = str2num(D{i,1});
    wormNc{i} = repmat(a,wormN(i,1),1);
end
wormNC = cell2mat(wormNc);
D = cell2mat(D(:,2));
D = [wormNC D];
D = D(:,ivar);
colname = regexprep(var,'[:]','_');
D = array2table(D,'VariableNames',colname);