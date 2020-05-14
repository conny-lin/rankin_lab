function StatT = ephys_randomSample(DataG,pSave,varargin)


%% setting
wormN = 60;
plateN = 100;
expN = 3;
vararginProcessor

%% process random samples ------------
%% random sampling

gN= size(DataG,2);


% decide group sequence
gu = cell(size(DataG,2),1);
for gi =1:size(DataG,2)
    gu{gi} = DataG(gi).name;
end
groupnameSet = gu([find(regexpcellout(gu,'N2'));find(~regexpcellout(gu,'N2'))]);
 
pValue = nan(expN,gN);
for ei = 1:expN
peakValue = nan(plateN,gN);
groupname = cell(size(peakValue));
peakTime = nan(size(peakValue));

for gi = 1:gN
    % get data
    D = DataG(gi).speedbm;
    time = DataG(gi).time(1,:);
    sampleN = size(D,1);
    gn = {DataG(gi).name};
    [~, gposition] = ismember(gn,groupnameSet);

for platei = 1:plateN    
    % baseline time
    tbs = (1:find(time==0)-1);
    tpeak = (find(time(1,:)==0)+1):find(time==2);
    % get random number
    DS = D(randi(sampleN,wormN,1),:);
    
    % find baseline (2SD above mean baseline)
    bsd = DS(:,tbs);
    [m,bs,~,sd,~] = ephys_findbaseline(bsd);
    
    % calculate peak
    t = time(tpeak);
    d = nanmean(DS(:,tpeak));
    [peakValue(platei,gposition), peakTime(platei,gposition)] = ephys_findrisepeak(d,t,bs);
    
    % record group name
    groupname(platei,gposition) = gn;
end

end

% stats
pValue(ei,:) = sum(~isnan(peakValue)) ./ repmat(plateN,1,gN);

end

%% flip for stats
g = repmat(groupnameSet',size(pValue,1),1);
g = reshape(g,numel(g),1);
x = reshape(pValue,numel(pValue),1);

%% 
StatT = grpstatsTable(x,g);
cd(pSave);
writetable(StatT,'ephys t28_30.csv');
%% sep group name
a = regexpcellout(g,'_','split');
strain = a(:,1);
dose = a(:,2);
dose(cellfun(@isempty,dose)) = {'0mM'};
gvarnames = {'strain','dose'}';
anovan_std(x,{strain dose},gvarnames,pSave,'ephys t28_30')


