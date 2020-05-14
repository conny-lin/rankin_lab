function [figure1, fignamePrefix] = PSTH(Data,varargin)

%% ADD FUNCTION PATH
pFun = {'/Users/connylin/Dropbox/Code/Matlab/Library/General';
        '/Users/connylin/Dropbox/Code/Matlab/Library/Graphs'};
for x = 1:numel(pFun); 
    p = pFun{x};
    if isdir(p) == 0; error('function path %d not found'); end
    addpath(p); 
end
%% setting
groupname = {};
mwtname =  {};
graphtitle = 'off';
summarytype = 'histogram';
visible = 0;
%% varargin processor
vararginProcessor;
%% generate recommended prefix
fignamePrefix = ['PSTH ',summarytype];
%% align time by first tap time
i = Data.tap == 1;
if sum(i) == 0; warning('no tap time; program %s abort',mfilename), figure1 = {}; return; end
taptime = [Data.id(i) Data.time(i)];
taptime = sortrows(taptime,1);
taptimeu = unique(Data.time(i));
i = [1;diff(taptime(:,1))]~=0;
tp = unique(taptime(i,2));
if numel(tp)~=1; error('tap time not standard'); end
tpstd = unique(floor(tp/10).*10);
Data.timestd = Data.time + (tpstd-tp);
% calculate tap line position
tpfloor = floor(unique(Data.time(Data.tap == 1)));
[~,linex] = ismember(taptimeu,unique(Data.time));
%% transform worm id into y axis
widu = unique(Data.id);
[~,j] = ismember(Data.id,widu);
Data.plot_wormid = j;
nworm = numel(widu);
%% make x axis (time)
timeu = unique(Data.timestd);
[~,j] = ismember(Data.timestd,timeu);
Data.plot_timeid = j;
ntime = numel(timeu);

%% make bias matrix
IMG = nan(nworm,ntime);
D = Data.bias;
for wi = 1:nworm
    i = Data.plot_wormid == wi;
    d = D(i);
    yi = Data.plot_timeid(i);
    IMG(wi,yi) = d;
end
Bias = IMG;
%% make speed matrix
IMG = nan(nworm,ntime);
D = Data.speed;
for wi = 1:nworm
    i = Data.plot_wormid == wi;
    d = D(i);
    yi = Data.plot_timeid(i);
    IMG(wi,yi) = d;
end
Speed = IMG;
%% make speed forward matrix
a = nan(size(Speed));
a(Bias== 1) = Speed(Bias==1);
% calculate forward average
b = a; b(isnan(b)) = 0;
SpeedFy = nanmean(b);        
a(isnan(Speed)) = Inf;
% create speed RGB map (red)
map = colormap_custom('grad1',0,'gradmax',0.9,'step',40,'color','r');
% scale to color map
b = colormap_scaledata(map,a);
SpeedFIMG = ind2rgb(b,map);
%% make speed rev matrix
a = nan(size(Speed));
a(Bias== -1) = Speed(Bias==-1);
b = a; b(isnan(b)) = 0;
SpeedRy = nanmean(b);
a(isnan(Speed)) = Inf;
% create speed RGB map (blue)
map = colormap_custom('grad1',0,'gradmax',0.9,'step',40,'color','b');
% scale to color map
 b = colormap_scaledata(map,a);
SpeedRIMG = ind2rgb(b,map);
%% make speed bias=0 matrix
a = nan(size(Speed));
a(Bias== 0) = Speed(Bias==0);
b = a; b(isnan(b)) = 0;
Speed0y = nanmean(b);
a(isnan(Speed)) = Inf;
% create speed RGB map (green)
map = colormap_custom('grad1',0,'gradmax',0.9,'step',40,'color','g');
% scale to color map
b = colormap_scaledata(map,a);
Speed0IMG = ind2rgb(b,map);
%% make speed bias=N/A matrix
a = nan(size(Speed));
a(isnan(Bias)) = Speed(isnan(Bias));
b = a; b(isnan(b)) = 0;
SpeedNaNy = nanmean(b);
a(isnan(Speed)) = Inf;
% create speed RGB map (grey)
map = colormap_custom('grad1',0,'gradmax',0.9,'step',40,'color','k');
% scale to color map
b = colormap_scaledata(map,a);
SpeedNaNIMG = ind2rgb(b,map);
%% make combine image
% spacer = ones(3,ntime,3);
IMG = ([SpeedRIMG;SpeedFIMG;Speed0IMG;SpeedNaNIMG]);
%% calculate summary plot
switch summarytype
    case 'mean speed'
        S = nan(size(Speed));
        S(Bias == -1) = -Speed(Bias==-1);
        S(Bias == 1) = Speed(Bias==1);
        S(Bias==0) = 0;
        Speedy = nanmean(S);
        PlotMatrix = [SpeedNaNy;Speed0y;SpeedFy; SpeedRy ;Speedy];
        plotcolor = {[.5 .5 .5],'g','r','b',[0.5 0 0.6]};
        ylabelname = 'Mean Speed';
        plotlegend = {'NA','Pause','Forward','Reverse','F+R'};
        ylimit = [min(min(PlotMatrix))*0.9 max(max(PlotMatrix))*1.05];
    case 'histogram'
        f = sum(Bias==1);
        r = sum(Bias==-1);
        p = sum(Bias==0);
        mia = sum(isnan(Bias));
        PlotMatrix = [mia;p;f; r];
        plotlegend = {'NA','Pause','Forward','Reverse'};
        plotcolor = {'k','g','r','b'};
        ylabelname = 'Bias Histogram';
        ylimit = [0 max(max(PlotMatrix))+1];
    case 'move speed'
        % calculate speed per category
        a = nan(size(Speed));
        a(isnan(Bias)) = Speed(isnan(Bias));
        SNaN = nanmean(a);
        a = nan(size(Speed));
        a(Bias==0) = Speed(Bias==0);
        S0 = nanmean(a);
        a = nan(size(Speed));
        a(Bias==1) = Speed(Bias==1);
        SF = nanmean(a);
        a = nan(size(Speed));
        a(Bias==-1) = Speed(Bias==-1);
        SR = nanmean(a);
        % calculate mean speed
        S = nan(size(Speed));
        S(Bias==-1) = -Speed(Bias==-1);
        S(Bias==1) = Speed(Bias==1);
        S(Bias==0) = 0;
        Speedy = nanmean(S);
        % create plotting matrix and graphing info
        PlotMatrix = [SNaN;S0;SF; SR ;Speedy];
        plotcolor = {[.5 .5 .5],'g','r','b',[0.5 0 0.6]};
        ylabelname = 'Mean Move Speed';
        plotlegend = {'NA','Pause','Forward','Reverse','F+R'};
        ylimit = [min(min(PlotMatrix))*0.9 max(max(PlotMatrix))*1.05];
end
%% PLOT
close;
if visible
    figure1 = figure('Color','none','Visible','on');
else
    figure1 = figure('Color','none','Visible','off'); 
end
hold on;
% curve plot on top
subplot1 = subplot(2,1,1,'Parent',figure1,'Box','off',...
    'Color','none','XTick',[],'XLim',[1 size(IMG,2)],'YLim',ylimit);
hold(subplot1,'all');
if isempty(groupname)==0 && isempty(mwtname)==0 && strcmp(graphtitle,'on')==1
    title(sprintf('%s [%s]',regexprep(groupname,'_',' '), regexprep(mwtname,'_','-')));
end
plot1 = plot(1:size(IMG,2),PlotMatrix,'Parent',subplot1,'LineWidth',2);
for si = 1:size(PlotMatrix,1)
    set(plot1(si),'Color',plotcolor{si},'DisplayName',plotlegend{si});
end
for xi = 1:numel(linex) % tap line
    x = linex(xi);
    line([x x],[0 nworm],'Color','k','LineWidth',1,'LineStyle',':')
end
ylabel(ylabelname);

% plot raster
axes1 = axes('Parent',figure1,'Position',[0.13 0.11 0.775 0.446552962298025],...
    'YDir','reverse','Box','off','Color','none','XLim',[1 size(IMG,2)],...
    'YLim',[1 size(IMG,1)],'YTick',[],'XTick',linex,'XTickLabel',tpfloor);
hold(axes1,'all');
image(IMG)
for xi = 1:numel(linex) % tap line
    x = linex(xi);
    line([x x],[0 size(IMG,1)],'Color','k','LineWidth',1,'LineStyle',':')
end
        
end























