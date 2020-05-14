function MWTSet = Dance_DrunkMoves_v1707(pMWT,pSave,varargin)
%% DRUNKMOVES
% inputs
%     pSave = '/Users/connylin/Dropbox/RL/Dance Output';
%     timeStartSet = 0;
%     timeIntSet = 60;
%     timeAssaySet = 60;
%     timeEndSet = Inf;


%% DEFAULTS & VARARGIN
outputsuffix = '';
timeStartSet = 0;
timeIntSet = 60;
timeAssaySet = 60;
timeEndSet = Inf;
saveimport = 0;
minSample = 15;
saveimport = 0;
vararginProcessor

%% STANDARD DANCE PROCESSING
% create MWTSet
[MWTSet,pSave] = DanceM_MWTSetStd_v1707(pMWT,mfilename('fullpath'),varargin,pSave);

%% specific processing
MWTSet.Info.timeAnalysisSet = [timeStartSet timeIntSet timeAssaySet timeEndSet];


%% CHOR
[MWTSet,MWTDB] = DanceM_chor(MWTSet,'DrunkPosture');


%% drunkposture2.dat
% IMPORT RAW DATA & COMPILE TO .MAT (.drunkposture2.dat)
Vind = MWTSet.Info.VarIndex;
[Import,MWTSet.Info.MsrList.drunkposture2] = import_drunkposture2(pMWT,MWTSet.MWTDB,MWTSet.Limit);

% validation
[Import,Import_bad] = validate_drunkposture2_endTime(Import,pMWT,Vind);
% store results
if saveimport == 1
    MWTSet.Import.drunkposture2 = Import;
    MWTSet.Import.drunkposture2_badendtime = Import_bad;
end


%% DETERMINE TIME ANALYSIS INTERVALS
fprintf('\nDetermining analysis times:\n');
D = Import;
% evaluate data time intervals
% get max start time
tStartThreshold = max(D.time(D.recordseq == 1));
mwtind = unique(D.mwtname)';
a = nan(size(mwtind));
for mwti = mwtind
    a(mwti) = max(D.time(D.mwtname == mwti));
end
tEndThreshold = min(a);
fprintf('-time threshold (start-finish): %.0f-%.0f(s)\n',...
    tStartThreshold,tEndThreshold);
% built assay times
if tStartThreshold > timeStartSet
    t1 = tStartThreshold;
else
    t1 = timeStartSet;
end
if timeEndSet > tEndThreshold; 
    tf = tEndThreshold; 
else
    tf = timeEndSet;
end
timeAnalysis = [t1 [t1+timeIntSet:timeIntSet:tf-timeIntSet]];
timeAnalysis(2,:) = (timeAnalysis+timeAssaySet);
timeAnalysis(2,1) = floor(timeAnalysis(2,1));
% add last block of time if remaining time is larger than half of the interval
if (tf - timeAnalysis(end,end)) > timeIntSet/2
    timeAnalysis = [timeAnalysis [timeAnalysis(end,end); tf]];
end
fprintf('-analysis blocks: %d\n',size(timeAnalysis,2));
MWTSet.Info.timeAnalysis = timeAnalysis;


%% CAL: MEAN OF EACH MSR PER TIME ANALYSIS INTERVAL PER PLATE
% get required variables
Msr = MWTSet.Info.MsrList.drunkposture2(~ismember(MWTSet.Info.MsrList.drunkposture2,'time'));
D = Import;
gnamesRef = MWTSet.Info.VarIndex.groupname;
mwtnameRef = MWTSet.Info.VarIndex.mwtname;
DbInd = MWTSet.Info.MWTDbInd;
timeAnalysis = MWTSet.Info.timeAnalysis;
% index time
D.timeind = nan(size(D,1),1);
for ti = 1:size(timeAnalysis,2)
    t1 = timeAnalysis(1,ti);
    t2 = timeAnalysis(2,ti);
    D.timeind(D.time >=t1 & D.time < t2) = ti;
end
% take out time not under analysis
D(isnan(D.timeind),:) = [];
% calculate
fprintf('\nCalculating descriptive stats per plate, this will take a while... ');
S = struct;
mwtind = unique(D.mwtname);
for mi = 1:numel(Msr)
    msr = Msr{mi};
    TS = table;
    for mwti = 1:numel(mwtind)
        info = repmat(DbInd(DbInd.mwtname == mwtind(mwti),:),size(timeAnalysis,2),1);
        d = D(D.mwtname == mwtind(mwti),:);
        [n,mn,sd,se,gn] = grpstats(d.(msr),d.timeind,{'numel','mean','std','sem','gname'});
        T = table;
        T.timeind = cellfun(@str2num,gn);
        T.assaytimeStart = timeAnalysis(1,:)';
        T.assaytimeFinish = timeAnalysis(2,:)';
        T.frame_N = n;
        T.mean = mn;
        T.SD = sd;
        T.SE = se;
        T = [info T];
        TS = [TS;T]; % store to central table;
    end
    S.(msr) = TS;
end
fprintf('Done\n');

% store
MWTSet.Data_Plate = S;


%% CAL: MSR MEAN PER GROUP (N=PLATES)
S = MWTSet.Data_Plate;
Msr = fieldnames(S);
A = struct;
for msri = 1:numel(Msr) 
    msr = Msr{msri};
    D = S.(msr);
    groupnameU = unique(D.groupname)';
    TS = table;
    for gi = groupnameU
        d = D(D.groupname == gi,:);
        [n,mn,sd,se,gn] = grpstats(d.mean, d.timeind,{'numel','mean','std','sem','gname'});
        T = table;
        T.groupname = repmat(gi,size(gn));
        T.timeind = cellfun(@str2num,gn);
        T.assaytimeStart = timeAnalysis(1,:)';
        T.assaytimeFinish = timeAnalysis(2,:)';
        T.N = n;
        T.mean = mn;
        T.SD = sd;
        T.SE = se;
        TS = [TS;T]; % store to central table;
    end
    A.(msr) = TS;
end
MWTSet.Data_GroupByPlate = A;


%% TRANSFORM DATA FOR GRAPH: CURVE
% prepare data for graphing
gnameR = MWTSet.Info.VarIndex.groupname;
% reorg table into graph
S = MWTSet.Data_GroupByPlate;
Msr = fieldnames(S);
A = struct;
for msri = 1:numel(Msr) 
    msr = Msr{msri};
    D = S.(msr);
    gnameU = unique(D.groupname)';
    tN = unique(D.timeind);
    K = struct;
    K.groupname = gnameR(gnameU);
    K.N = nan(size(gnameU));
    B = nan(numel(tN),numel(gnameU));
    K.timeind = B;
    K.mean = B;
    K.SE = B;
    for gi = 1:numel(gnameU)
        d = D(D.groupname == gnameU(gi),:);
        K.N(gi) = mean(d.N);
        K.timeind(d.timeind,gi) = d.timeind;
        K.mean(d.timeind,gi) = d.mean;
        K.SE(d.timeind,gi) = d.SE;
    end
    A.(msr) = K;
end
MWTSet.Graph.Curve = A;


%% GRAPH: CURVE - INDIVIDUAL
pSaveA = [pSave,'/Graph Curve']; if isdir(pSaveA) == 0; mkdir(pSaveA); end
G = MWTSet.Graph.Curve;
color = [0 0 0;1 0 0;[0.04 0.52 0.78];[0.478 0.0627 0.8941]]; 
Msr = fieldnames(G);
for mi = 1:numel(Msr)
    msr = Msr{mi};
    X = G.(msr).timeind;
    Y = G.(msr).mean;
    E = G.(msr).SE;
    gname = MWTSet.Graph.Curve.(msr).groupname;
    PlateN = MWTSet.Graph.Curve.(msr).N;
    titlestr = gen_Nstring(PlateN);
    figure1 = Graph_errorbar(X,Y,E,gname,color,...
        'titlestr',titlestr,'xname','time(min)','yname',msr,'visiblesetting','off');
    savefigepsOnly150(sprintf('%s',msr),pSaveA)
end


%% GRAPH: CURVE - COMPOSITE
fprintf('Generating composite graph\n');
pSaveA = pSave;
G = MWTSet.Graph.Curve;
Msr = {'number' 'goodnumber' 'speed' 'bias',...
    'midline' 'area' 'width' 'morphwidth',...
    'length' 'aspect' 'kink' 'curve'};
nM = numel(Msr);
nMs = sqrt(nM);
nGY = floor(nMs);
nGX = ceil(nMs);
gname = MWTSet.Graph.Curve.(msr).groupname;
color = [0 0 0;1 0 0;[0.04 0.52 0.78];[0.478 0.0627 0.8941]]; 
PlateN = MWTSet.Graph.Curve.(msr).N;
titlestr = gen_Nstring(PlateN);
xname = 'time(min)';
% make figure
figure1 = figure('Color',[1 1 1],'Visible','off');

for mi = 1:numel(Msr)
    msr = Msr{mi};
    X = G.(msr).timeind;
    Y = G.(msr).mean;
    E = G.(msr).SE;
    % Create subplot
    subplot1 = subplot(nGY,nGX,mi,'Parent',figure1,'XTick',0:10:max(max(X)),...
        'FontSize',8);
    xmax = max(max(X));
    xlim(subplot1,[0 xmax+1]);
    % ylim(subplot1,[0 150]);
    hold(subplot1,'all');
    errorbar1 = errorbar(X,Y,E);
    if size(color,1) > size(Y,2)
        colorn = size(Y,2);
    else
        colorn = size(color,1);
    end
    for x = 1:colorn
        set(errorbar1(x),'Color',color(x,1:3))
    end
    for x = 1:size(Y,2)
        set(errorbar1(x),'DisplayName',regexprep(gname{x},'_',' '))
    end
    title(msr);
    % put x lable at the bottom row only
    if mi > nGX*((nGY)-1)
        xlabel(xname,'FontSize',8); % Create xlabel
    end
    % ylabel(msr); % Create ylabel
    % legend
    if mi == 1
        legend1 = legend(subplot1,'show');
        set(legend1,'EdgeColor',[1 1 1],...
            'Position',[0.02 0.88 0.025 0.08],...
            'FontSize',6); 
    end
end
% Create textbox
annotation(figure1,'textbox',...
    [0.01 0.99 0.33 0.018],...
    'FitBoxToText','off','String',titlestr,...
    'LineStyle','none','FontSize',8);
savefigepsOnly150('Curve composit',pSaveA)


%% Export curve data
pSaveA = [pSave,'/Data Curve'];
export_drunkposture2_bygroup(MWTSet.Graph.Curve,pSaveA)



%% Export
cd(pSave);
save(mfilename,'MWTSet');




