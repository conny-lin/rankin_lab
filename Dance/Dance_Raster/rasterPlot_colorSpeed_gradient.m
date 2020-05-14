function figure1 = rasterPlot_colorSpeed_gradient(Data,time,visible)
%%

%% convert speed to graphing gradient
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


%% create gradient based on gradmax
SpeedGrad = nan(size(Data));

i = Data > 0;
F = Data./gradmax;
SpeedGrad(i) = F(i);

i = Data < 0;
R = Data./-gradmin;
SpeedGrad(i) = R(i);

i = Data == 0;
SpeedGrad(i) = 0;

SpeedGrad(SpeedGrad > 1) = 1;
SpeedGrad(SpeedGrad < -1) = -1;


%% create time label
% time = rTime(1,:);
timeticks = 1:10:numel(time)-10;
timelabel = time(timeticks);
timelabel = round(timelabel);
timelabel = regexprep(cellstr(num2str(timelabel')),' ','')';


%% plot per worm
Img = SpeedGrad;
if visible
    figure1 = figure('Color',[1 1 1],'Visible','on');
else
    figure1 = figure('Color',[1 1 1],'Visible','off');
end
axes1 = axes('Parent',figure1,...
    'YDir','reverse',...
    'YColor',[0.8 0.8 0.8],...
    'XColor',[0.8 0.8 0.8],...
    'XTickLabel',timelabel,...
    'XTick',timeticks,...
    'Layer','top');
imagesc(Img,'Parent',axes1)
