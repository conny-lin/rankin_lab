function fig_bar(X,Y,E,GN,genotype,pSave,varargin)
% 
% %% graphic setting --------
color = [0 0 0; 0 0 0; 1 0 0; 1 0 0];
colorface = {'none';[0 0 0];'none';[1 0 0]};
marker = {'o','o','o','o'};
linestyle = {'-','-','-','-'};
time = 0.1:0.1:0.5;
ylim = [0 1];
xlim = [0.5 4.5];
yname = '';
dosenames = {'0mM','400mM','0mM','400mM'};
ylable = '';
w = 3;
h = 3;
vararginProcessor;



%% plot bar graphs ----------------------
close;
figure1 = figure('Visible','off',...
                'PaperUnits','inches',...
                'PaperPosition',[0 0 w h]);
            
            
% Create axes
axes1 = axes('Parent',figure1,...
            'XTick',1:numel(X),...
            'XLim',xlim,...
            'XTickLabel',dosenames,...
            'YLim',ylim,...
            'Position',[0.15 0.15 0.85 0.84]);       
ylabel(yname);
hold(axes1,'on');
% ylabel(ylabel);
% create bar graph
for gi = 1:numel(GN)
    bar('Parent',axes1,X(gi),Y(gi),...
        'FaceColor',colorface{gi},...
        'EdgeColor',color(gi,:),...
        'LineStyle',linestyle{gi});
end


% Create multiple error bars using matrix input to errorbar
errorbar1 = errorbar('Parent',axes1,X,Y,E);
set(errorbar1,'LineStyle','none','Color',[0 0 0])

% return
% for gi = 1:numel(gn)
%     set(errorbar1(gi),...
%         'Color',color(gi,:),...
%         'Marker','none',...
%         'MarkerEdgeColor',color(gi,:),...
%         'MarkerFaceColor',colorface{gi},...
%         'LineStyle','none');
% end
% return

%  create text box and annotations
height_text = 0.08;
% Create textbox
annotation(figure1,'textbox',...
    [0.3 height_text 0.1 0],...
    'String',genotype{1},...
    'HorizontalAlignment','center',...
    'FitBoxToText','off',...
    'EdgeColor','none');

% Create textbox
annotation(figure1,'textbox',...
    [0.75 height_text 0.1 0],...
    'String',genotype{2},...
    'HorizontalAlignment','center',...
    'FitBoxToText','off',...
    'EdgeColor','none');

% 
% % Create line
% height_line= 0.07;
% annotation(figure1,'line',[0.18 0.45],[height_line height_line]);
% 
% % Create line
% annotation(figure1,'line',[0.61 0.88],[height_line height_line]);
% 





% save figure ------------------------------------
print(pSave,'-dtiff'); 
[pS,sname] = fileparts(pSave);
printfig(sname,pS,'closefig',1,'w',w,'h',h)
% -------------------------------------------------

























