%% OBJECTIVE:
% find genes can be analyzed in 10sISI and group into batches



%% paths & settings
pAnalysis = '/Users/connylin/Dropbox/Lab/MWT_Analysis';
% function path
addpath('/Users/connylin/Dropbox/MATLAB/Function_Library_Public');

%% survey data 
[~,~,Ef,pE] = dircontent(pAnalysis);
% find strains that is 10sISI
i = regexpcellout(Ef,'100s30x10s');
Ef = Ef(i); pE = pE(i);
[~,~,fG,pG] = cellfun(@dircontent,pExp,'UniformOutput',0);
pG = celltakeout(pG); fG = celltakeout(fG);
a = regexpcellout(fG,'_','split');
strain = a(:,1);
strainU = unique(strain);
% find strains have alcohol experiments
a = cellfun(@strcat,strainU,cellfunexpr(strainU,'_400mM'),'UniformOutput',0);
strainU = unique(strain(ismember(fG,a)));
a = cellfun(@strcat,strainU,cellfunexpr(strainU,'_400mM'),'UniformOutput',0);
groupU = [strainU; a]; % find all groups
char(strainU)
T = table;
T.strain = strainU;
cd('/Users/connylin/Downloads');
writetable(T,'strain.csv');


%% get npr-1 list
cd('/Users/connylin/Dropbox/Lab/MWT_Analysis_ByExp');
strainList = readtable('/Users/connylin/Dropbox/Lab/MWT_Analysis_ByExp/strain_pathway.csv');
npr1List = strainList(ismember(strainList.pathway,'npr-1'),:);






return


%% load data
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
    