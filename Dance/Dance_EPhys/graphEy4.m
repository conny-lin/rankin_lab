function graphEy4(G,genotype,strain,timeName,pSave,Setting,varargin)

%% defults
zeroline = false;
switch timeName
    case 't1'
        ylimNum = [-.5 .4];
    case 't28_30'
        ylimNum = [-.3 .4];
end 
xlimNum = [-0.15 1];
vararginProcessor

%% figure
figure1 = figure('Visible','on'); hold on;

%% zeroline
if zeroline
line(xlimNum,[0 0],'Color',[0.6 0.6 0.6],'LineStyle','-','Linewidth',1)
end

%% errorbar
errorbar1 = errorbar(G.x,G.y,G.e);

%% place setting
SettingNames = fieldnames(Setting);
for si = 1:numel(SettingNames)
    nms = SettingNames{si};
for gi = 1:size(errorbar1,2)
    errorbar1(gi).(nms) = Setting.(nms){gi};
end
end
xlim(xlimNum);
ylim(ylimNum);



%% add descriptive components
title(genotype,'FontSize',8)
ylabel('velocity (body length/s)','FontSize',9)



%% save
savename = sprintf('%s ephys %s',strain,timeName);
printfig(savename,pSave,'w',2.5,'h',2.5,'closefig',0,'version',3);



