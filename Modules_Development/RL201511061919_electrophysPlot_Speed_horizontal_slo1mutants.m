%% HYPOTHESIS: 
%     after tap, (at later tap) worms do forward > short reverse > pause
%     alcohol knocks out this resposne
%     find, for worms who does pause/forward > reverse > pause 
%     pause/forward > reverse (mean duration) > pause (mean duration)



%% paths & settings
pHome = '/Users/connylin/Dropbox/Lab/MWT_Analysis_ByExp/slo1 10sISI';
frameInt = 0.1;
NWORMS = Inf;
alphaValue = 0.05;
posttestname = 'bonferroni';
% time
tapTimeList = [100:10:100+10*29];
tpre = tapTimeList- 3;
tpost = tapTimeList + 5;
groupTarget = {'N2','N2_400mM','NM1968','NM1968_400mM'};

% exp info path
pSave = sprintf('%s/Stats/tapRShiftGraph',pHome);
if isdir(pSave) == 0; mkdir(pSave); end
pAnalysis = '/Users/connylin/Dropbox/Lab/MWT_Analysis';
% function path
addpath('/Users/connylin/Dropbox/MATLAB/Function_Library_Public');
addpath('/Users/connylin/Dropbox/MATLAB/programs_rankinlab/Library/Modules/rasterPlot_colorSpeed');


%% get groups
% load global data
plateInfo = readtable(sprintf('%s/plate_info.csv',pHome));
% narrow down group name target
plateInfo(~ismember(plateInfo.groupname,groupTarget),:) = [];
% create output folder
g = unique(plateInfo.groupname);
pG = cellfun(@strcat,cellfunexpr(g,pHome),cellfunexpr(g,'/'),g,'UniformOutput',0);


% includ only MWT have necessary files for tap times
A = plateInfo;
display 'Group plates found:';
tabulate(plateInfo.groupname)
display ' ';
val = true(size(plateInfo,1),1);
for x = 1:size(plateInfo,1)
    pMWT = sprintf('%s/%s/%s/%s',pAnalysis,A.expname{x},A.groupname{x}, A.mwtname{x});
    [~,pset] = dircontent(pMWT,'*.set');
    if isempty(pset) == 1
        val(x) = false;
    else
        p = char(pset);
        [~,p] = dircontent(pMWTp,[p(1:regexp(p,'(.set)')-1),'.dat']);
        if isempty(p) == 1
            val(x) = false;
        end
    end
end
plateInfoOld = A;
plateInfo(~val,:) = [];
% report on size
display 'Group plates with valid tap .dat files:';
tabulate(plateInfo.groupname)
display ' ';



%% calculation
gname = {};
% start
fprintf('calculating data...\n');
Data = {};
for ti = 1:numel(tapTimeList)
    start = tpre(ti); 
    finish = tpost(ti);
    taptime = tapTimeList(ti);
    tapCol = ((taptime-start)/frameInt)+1;
    frameNumberExpected = numel((start:frameInt:finish-frameInt));
    fprintf('time(%d/%d): %d-%ds\n',ti,numel(tapTimeList),start,finish);
    for gi = 1:numel(pG)
        %% get dir to MWT
        [~,fG] = fileparts(pG{gi});
        gname{gi}= fG;
        fprintf('-%s\n',fG);
        i = find(ismember(plateInfo.groupname,fG));
        pMWT = cell(numel(i),1);
        A = plateInfo(i,:);
        for x = 1:numel(i)
            p = sprintf('%s/%s/%s/%s',pAnalysis,A.expname{x},A.groupname{x}, A.mwtname{x});
            pMWT{x} = p;
        end

        %% get trinity data
        D = [];
        for mwti = 1:numel(pMWT)
            pMWTp = pMWT{mwti};
            addpath('/Users/connylin/Dropbox/MATLAB/programs_rankinlab/Library/Modules/rasterPlot_colorSpeed');
            a = trinitySummary(pMWTp,start,finish,'frameInt',frameInt,'NWORM',Inf,'saveoption',0);        
            [~,b] = fileparts(pMWTp);
            if isempty(a) == 0
                % if a is bigger than D, which happens when tap lands on a integer such as 100
                if size(a,2) > frameNumberExpected
                    a = a(:,1:frameNumberExpected);
                end
                if isempty(D) ==1; D = a; else D = [D;a]; end
            else
                fprintf('-- no valid worms for [%s] at time %d-%ds\n',b,start,finish);
                
            end
        end

        %% mean response rate
        if gi == 1
            M = mean(D);
            N = size(D,1);
            SE = std(D)./sqrt(size(D,1)-1);
        else
            M = [M;mean(D)];
            N = [N;size(D,1)];
            SE = [SE;std(D)./sqrt(size(D,1)-1)];
        end

    end
    fprintf('\n');

    % get time point summary
    Data{ti}.M = M;
    Data{ti}.N = N;
    Data{ti}.SE = SE;
end
cd(pSave); save('temp.mat','Data');

%% Get Mean speed per group
% find out how many groups
gN = size(pG);
% find out group names
[~,gnames] = cellfun(@fileparts,pG,'UniformOutput',0);
DataG = cell(numel(gnames),2);
rowN = size(Data{1}.M,2);
for g = 1:numel(gnames)
for x = 1:numel(Data)
    DataG{g,2} = gnames{g};
    DataG{g,1}(x,1:rowN) = Data{x}.M(g,1:rowN);
end
end

%% get y min/max (spacer of vertical axis) 
yGraphSpacerRatio = 1.01;
for g = 1:numel(gnames)
    a = min(min(DataG{g,1}));
    b = max(max(DataG{g,1}));
    if g == 1
        minY = a;
        maxY = b;
    else
        if minY > a; minY = a;end
        if maxY < b; maxY = b; end
    end
end


% calculate total space
if (minY < 0 && maxY < 0) || (minY > 0 && maxY > 0)
    ySize = maxY - minY;
elseif minY < 0 && maxY > 0
    ySize = abs(minY) + maxY;
else
    error('some possibilities not considered, code');
end
a = ySize*yGraphSpacerRatio;
yspacer = 0:a:numel(gnames)*a-a;
if numel(yspacer) ~= numel(gnames); 
    error('yspacer not equal to group'); 
end


%  x steps (spacer for horizontal axis)
xBlankSize = 10;
xgraphN = size(DataG{1,1},1);
nX = size(DataG{1,1},2);
xGraphSpacer = nX+xBlankSize;
xstand = 1:nX;
xspacer = 0:xGraphSpacer:xgraphN*xGraphSpacer-xGraphSpacer;
if numel(xspacer) ~= xgraphN; 
    error('xspacer not equal to data point number'); 
end
xMatrix = nan(size(DataG{1,1}));
for i = 1:xgraphN
    xMatrix(i,:) = xstand + xspacer(i);
end
if numel(unique(xMatrix)) ~= numel(xMatrix)
    error('xspacer overlap');
end


% plot
% plot setting
gColor = [0 0 0; 1 0 0; 0.87058824300766 0.490196079015732 0;...
    1 0.800000011920929 0.400000005960464;...
    0 1 0; 0 0 1; 0.47843137383461 0.062745101749897 0.894117653369904];
% create figure
figure1 = figure('Color',[1 1 1]);
axes1 = axes('Parent',figure1,'ZColor',[1 1 1],'YColor',[1 1 1],...
    'XColor',[1 1 1]);
box(axes1,'on');
hold(axes1,'all');
for gi = 1:numel(gnames)
    d = DataG{gi,1} - yspacer(gi);
    for ti = 1:size(d,1)
        plot(xMatrix(ti,:),d(ti,:),'Color',gColor(gi,1:3),'LineWidth',1)
        hold on
    end
end
% draw zero lines (=yspacer)
xmax = max(xspacer)+max(xstand);
for x = 1:numel(yspacer)
    line([0 xmax], [-yspacer(x) -yspacer(x)],'Color',[0 0 0])
end
xlim([0 max(max(xMatrix))+1]);
ylim([min(min(d)) 0+maxY]);

% save
h = (gcf);
set(h,'PaperSize',[xgraphN numel(gnames)]); % set to save as appeared on screen
set(h,'PaperPosition',[.1 .1 xgraphN numel(gnames)+.001]);
titlename = sprintf('%s/electrophys_speed',pSave);
print (h,'-dpdf', '-r0', titlename); % save as (r600 will give better resolution)
close all;


return
    