function [figure1,Img] = rasterPlot_colorSpeed_gradient2(Data,time,visible,varargin)

%% defaults
axeslabel = false;
vararginProcessor

%% define max color gradient
% define max min
gradmax = 0.6;
gradmin = -0.8;
% check max min 3SD speed
d = reshape(Data,numel(Data),1);
f = d(d > 0);
speed_max = mean(f) + 3*std(f);
if speed_max > gradmax; warning('forward+3SD speed > %.1f',gradmax); end
r = d(d < 0);
speed_min = mean(r) - 3*std(r);
if speed_min < gradmin; warning('reverse-3SD speed < %.1f',gradmin); end


%% convert speed 2 color based on gradmax
Speedcolor = nan(size(Data));

i = Data > 0;
F = Data./gradmax;
Speedcolor(i) = F(i);

i = Data < 0;
R = Data./-gradmin;
Speedcolor(i) = R(i);

i = Data == 0;
Speedcolor(i) = 0;

Speedcolor(Speedcolor > 1) = 1;
Speedcolor(Speedcolor < -1) = -1;


%% create time label
% time = rTime(1,:);
timeticks = 1:10:numel(time)-10;
timelabel = time(timeticks);
timelabel = round(timelabel);
timelabel = regexprep(cellstr(num2str(timelabel')),' ','')';


%% plot per worm
Img = Speedcolor;
if visible
    figure1 = figure('Color',[1 1 1],'Visible','on');
else
    figure1 = figure('Color',[1 1 1],'Visible','off');
end
% set paper size
set(figure1, 'PaperPosition',[0 0 2 3])

colormap('jet');
axes1 = axes('Parent',figure1);
hold(axes1,'on');

% set axes
% axes1 = axes('Parent',figure1,...
%     'YColor',[1 1 1],...
%     'XColor',[1 1 1],...
%     'XTickLabel',timelabel,...
%     'XTick',timeticks,...
%     'Layer','top');

% expand graph to fill the page
sz = size(Speedcolor);
xlim([0 sz(2)])
ylim([0 sz(1)])


% remove axes
if ~axeslabel
    set(axes1,'Layer','top','XColor','none','YColor','none');
end

imagesc(Img,'Parent',axes1)






































