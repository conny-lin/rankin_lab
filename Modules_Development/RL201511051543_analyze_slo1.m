%% OBJECTIVE:
% find genes can be analyzed in 10sISI and group into batches



%% paths & settings
frameInt = 0.1;
NWORMS = Inf;
alphaValue = 0.05;
posttestname = 'bonferroni';
% time
tapTimeList = [100:10:100+10*29];
tpre = tapTimeList- 1;
tpost = tapTimeList + 5;

pAnalysis = '/Users/connylin/Dropbox/Lab/MWT_Analysis';
pHome = '/Users/connylin/Dropbox/Lab/MWT_Analysis_ByExp/slo1 10sISI';
% function path
addpath('/Users/connylin/Dropbox/MATLAB/Function_Library_Public');
pSave = sprintf('%s/Stats/tapRShiftGraph',pHome);
if isdir(pSave) == 0; mkdir(pSave); end



%% get list and paths
cd('/Users/connylin/Dropbox/Lab/MWT_Analysis_ByExp');
strainList = readtable('/Users/connylin/Dropbox/Lab/MWT_Analysis_ByExp/strain_pathway.csv');
strainList = strainList(ismember(strainList.pathway,'slo-1'),:);
% get 400mM
s = strainList.strain;
groupT = [s;cellfun(@strcat,s,cellfunexpr(s,'_400mM'),'UniformOutput',0)];
% load global data
plateInfo = readtable(sprintf('%s/plate_info.csv',pHome));
% get only info in strain list
g = plateInfo.groupname;
plateInfo = plateInfo(ismember(g,groupT),:);



%% get path to home analysis group folder
pG = cellfun(@strcat,cellfunexpr(groupT,pHome), cellfunexpr(groupT,'/'),groupT,'UniformOutput',0);
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
    fprintf('time(%d/%d): %d-%ds',ti,numel(tapTimeList),start,finish);
    for gi = 1:numel(pG)
        %% get dir to MWT
        [~,fG] = fileparts(pG{gi});
        gname{gi}= fG;
        fprintf(', %s',fG);
        i = find(ismember(plateInfo.groupname,fG));
        pMWT = cell(numel(i),1);
        A = plateInfo(i,:);
        for x = 1:numel(i)
            p = sprintf('%s/%s/%s/%s',pAnalysis,A.expname{x},A.groupname{x}, A.mwtname{x});
            pMWT{x} = p;
        end

        %% get trinity data
        for mwti = 1:numel(pMWT)
            pMWTp = pMWT{mwti};
            addpath('/Users/connylin/Dropbox/MATLAB/programs_rankinlab/Library/Modules/rasterPlot_colorSpeed');
            a = trinitySummary(pMWTp,start,finish,'frameInt',frameInt,'NWORM',Inf,'saveoption',0);        
            % if a is bigger than D, which happens when tap lands on a integer such as 100
            if size(a,2) > frameNumberExpected
                a = a(:,1:frameNumberExpected);
            end
            if mwti == 1; D = a; else D = [D;a]; end
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




%% electrophys plot load data
D = load('/Users/connylin/Dropbox/Lab/MWT_Analysis_ByExp/Dose_10sISI/Stats/tapRShiftGraph/temp.mat');
Data = D.Data;

% reorganize
M = nan(numel(Data),60);
for ti = 1:numel(Data)
    M(ti,:) = Data{ti}.M(fGi,:);
end

% x spacer
xstand = [1:size(M,2)];
rowNgraph = 5;
xspacer = 0:max(xstand)+10:max(xstand)*rowNgraph+rowNgraph;

% y spacer
a = nan(numel(Data)/5,1);
b = a;
n = 0;
list = 1:5:size(M,1)-4;
for ti = list
    n = n+1;
    d = M(ti:ti+4,:);
    a(n) = min(min(d));
    b(n) = max(max(d));
end
if max(a) < 0
    yint = (max(b)-max(a))*1.5;
elseif max(a) > 0
    yint = (max(b))*2;
end
c = [0:yint:yint*((numel(Data)/5-1))];
yspacer = cumsum([0;a(2:end)]) - c(1:numel(list))';

%% yspacer ***** pre defined by N2_100mM data ******
yspacer = [0   -0.6995   -1.4016   -2.0277   -2.6432   -3.2089];


%% create figure
% plot setting
gColor = [0 0 0; 1 0 0; 0.87058824300766 0.490196079015732 0;...
    1 0.800000011920929 0.400000005960464;...
    0 1 0; 0 0 1; 0.47843137383461 0.062745101749897 0.894117653369904];

figure1 = figure('Color',[1 1 1]);
axes1 = axes('Parent',figure1,'ZColor',[1 1 1],'YColor',[1 1 1],...
    'XColor',[1 1 1]);
box(axes1,'on');
hold(axes1,'all');
for ti = 1:numel(list)
    d = M(ti:ti+4,:);
    for gi = 1:size(d,1)
        x = xstand + xspacer(gi);
        y = d(gi,:) + yspacer(ti);
        plot(x,y,'Color',gColor(1,1:3),'LineWidth',2)
        hold on
    end
end

% draw zero lines
xmax = max(xspacer)+max(xstand);
for x = 1:numel(yspacer)
    line([0 xmax], [yspacer(x) yspacer(x)],'Color',[0.5 0.5 0.5])
end

title({regexprep(fG,'_',' ')},'FontSize',16);

% save
h = (gcf);
set(h,'PaperSize',[8.5 11]); % set to save as appeared on screen
set(h,'PaperPosition',[.25 .25 [8.5 11]-0.25]);
titlename = sprintf('%s/%s',pSave,fG);
print (h,'-dpdf', '-r0', titlename); % save as (r600 will give better resolution)
close all;


return
    