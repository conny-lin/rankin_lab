function [fig,ax1,e1] = graphHabCurvePck(X,Y,E,msr,gnu,varargin)

% varargin setting ++++++++++++++
% figure size
w = 4;
h = 3;
% visibility
visibleopt = 'off';
genotype = '';
graphPack = 'cathyline';
markersize = 4.5;
vararginProcessor
% --------------------------

% get graph pack +++++++++++++
gp = graphsetpack(graphPack);
gpn = fieldnames(gp);
% ----------------------------

% generate graphic variables +++++++++
% y axis settings
switch msr
    case 'RevFreq'
        yaxislabel = 'P (reversal)';
        ylimarray = [0 1];
    otherwise
        yaxislabel = msr;
        y1 = Y+E;
        ymax = max(max(y1))*1.05;
        ymax = round(ymax,2);
        y2 = Y-E;
        ymin = min(min(y2));
        ymin = ymin - ymin*0.05;
        ymin = round(ymin,2);
        ylimarray = [ymin ymax];
end
% x axis settings
xtickrange = 0:10:30;
xlimarray = [0 30.5];
xaxislabel = 'Tap';
% prepare gns names
gnss = regexprep(gnu,'_',' ');
% ----------------------------------

% plot structure ++++++++++++++++
% figure
fig = figure('Visible',visibleopt,'PaperSize',[w h],...
            'Unit','Inches',...
            'Position',[0 0 w h]); 
% axes
ax1 = axes('Parent',fig,'XTick',xtickrange); 
hold(ax1,'on');
% axis setting
ax1.TickLength = [0 0];
ax1.FontSize = 10;
% title
title(genotype)
% x axis 
xlabel(xaxislabel)
xlim(xlimarray);
% y axis 
ylim(ylimarray);
ylabel(yaxislabel);
% -------------------------------

% plot variables +++++++++++++++
e1 = errorbar(X,Y,E,'Marker','o','MarkerSize',markersize);
% apply settings 
for gi = 1:numel(gp.Color)
    for gpi = 1:numel(gpn)
        e1(gi).(gpn{gpi}) = gp.(gpn{gpi}){gi};
    end
end
for gi = 1:numel(gnu)
   e1(gi).DisplayName = gnu{gi};
end
% -----------------------------

% legend ------------------------
%     lg = legend('show');
%     lg.Location = 'eastoutside';
%     lg.Box = 'off';
%     lg.Color = 'none';
%     lg.EdgeColor = 'none';
%     lg.FontSize = 8;
% -------------------------------












