function DataG = ephys_tranformData2struct(Data,gu)

%%
gidu = unique(Data.groupid);
DataG = struct;
for gi = 1:numel(gidu)
%     processIntervalReporter(numel(gidu),1,'group:',gi); % interval reporter
    Data1 = Data(Data.groupid==gidu(gi),:); % get data
    DataG(gi).name = gu{gi}; % record name
   
    % separate by worms
   tu = unique(Data1.timeround); % unique time
   ID = unique([Data1.mwtid Data1.timeid Data1.id],'rows'); % create id combo
   S = nan(size(ID,1),numel(tu)); % declare speed array
   B = S; % declare speed/midline array
   T = S; % declare time array
   % create index matrix
   matrixSize = size(T);
   sourceInd = nan(size(Data1,1),1);
   % summarize per time point
   for ti = 1:numel(tu)
       processIntervalReporter(numel(tu),20,' time:',ti)
       i = Data1.timeround == tu(ti);
       D = Data1(i,:);
       [j,k] = ismember(ID,[D.mwtid D.timeid D.id],'rows');
       linearInd = sub2ind(matrixSize, k(j),repmat(ti,sum(j),1));
       if numel(linearInd) ~= sum(i)
           error('bad');
       end
       sourceInd(i)=linearInd;
   end
   clear D;
   T(sourceInd) = Data1.timeround;
   S(sourceInd) = Data1.speed.*Data1.bias;
   B(sourceInd) = (Data1.speed.*Data1.bias)./mean(Data1.midline);
   
   
   % store in structural array
   DataG(gi).speedb = S;
   DataG(gi).speedbm = B;
   DataG(gi).time = T;
   % add id (20161002v)
   DataG(gi).id = array2table(ID,'VariableNames',{'mwtid','timeid','wormid'});
end