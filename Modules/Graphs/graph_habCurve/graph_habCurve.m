function graph_habCurve(X,Y,E,GroupName,pSaveA,yaxislabel,savename,plateN,titlename,color)

if nargin < 10
    color = [0 0 0; 1 0 0; [0.5 0.5 0.5]; [0.04 0.52 0.78]]; 
end

% Create figure
figure1 = figure('Color',[1 1 1],'Visible','off');
axes1 = axes('Parent',figure1,'FontSize',18);
box(axes1,'off');
xlim(axes1,[0 size(Y,1)+1]);
ylim(axes1,[min(min(Y-E))*.9 max(max(Y+E))*1.1]);
hold(axes1,'all');
errorbar1 = errorbar(X,Y,E,'LineWidth',1);
n = min([size(color,1) size(Y,2)]);
for x = 1:n
    set(errorbar1(x),'Color',color(x,1:3))
end

for x = 1:size(Y,2)
    set(errorbar1(x),'DisplayName',GroupName{x})
end
title(titlename,'FontSize',12);
xlabel('Tap','FontSize',18);
ylabel(yaxislabel,'FontSize',18);
legend1 = legend(axes1,'show');
set(legend1,'EdgeColor',[1 1 1],...
    'Location','NorthEastOutside',...
    'YColor',[1 1 1],'XColor',[1 1 1],'FontSize',12);

% text strings
a = plateN;
b = '';
for x = 1:numel(a)
    b = [b,num2str(a(x)),', '];
end
str1 = ['N(plate) = ',b(:,1:end-2)];
textboxstr = [str1];
annotation(figure1,'textbox',...
    [0.70 0.015 0.256 0.05],...
    'String',textboxstr,...
    'FitBoxToText','off',...
    'EdgeColor',[1 1 1]);
% save fig
savefigepsOnly(savename,pSaveA);



