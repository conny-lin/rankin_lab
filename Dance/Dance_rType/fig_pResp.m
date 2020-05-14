function fig_pResp(pM,X,Y,E,titlename,gn,varargin)

dispname = {'Acceleration','Reversal','Pauses'};
vararginProcessor

% figure ---------------------------------------
% Create figure
figure1 = figure('Visible','off','Position',[0 0 3 4]);

% Create axes
axes1 = axes('Parent',figure1);
hold(axes1,'on');

% Create multiple error bars using matrix input to errorbar
errorbar1 = errorbar(X,Y,E,'Marker','o');
color = [[1 0 0]; [0 0.447058826684952 0.74117648601532]; [0 0 0]];
for ei = 1:size(X,2)
set(errorbar1(ei),'DisplayName',dispname{ei},...
        'MarkerFaceColor',color(ei,:),...
    'MarkerEdgeColor',color(ei,:),...
    'Color',color(ei,:));

end

% Create xlabel
xlabel('time(s)');
% Create title
title(titlename);
% Create ylabel
ylabel('P (responses)');
% Uncomment the following line to preserve the X-limits of the axes
xlim(axes1,[X(1) X(end)]);
% Uncomment the following line to preserve the Y-limits of the axes
ylim(axes1,[0 1]);
% Create legend
legend1 = legend(axes1,'show');
set(legend1,'Location','north','EdgeColor',[1 1 1]);

% save
savename = [gn,' no pauses'];
printfig(savename,pM,'w',5,'h',4)
end