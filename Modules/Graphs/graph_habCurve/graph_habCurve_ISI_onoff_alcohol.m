function graph_habCurve_ISI_onoff_alcohol(YMatrix1, EMatrix1,msrname)
%CREATEFIGURE(YMATRIX1, EMATRIX1)
%  YMATRIX1:  errorbar y matrix
%  EMATRIX1:  errorbar e matrix

%  Auto-generated by MATLAB on 31-Jan-2016 15:28:36

% Create figure
figure1 = figure('Color',[1 1 1]);

% Create axes
axes1 = axes('Parent',figure1,'FontSize',24);
%% Uncomment the following line to preserve the X-limits of the axes
xlim(axes1,[0 31]);
hold(axes1,'all');

% Create multiple error bars using matrix input to errorbar
errorbar1 = errorbar(YMatrix1,EMatrix1,'MarkerSize',5,'Marker','o',...
    'LineWidth',2);
set(errorbar1(2),'MarkerFaceColor',[1 0 0],'MarkerEdgeColor',[1 0 0],...
    'DisplayName','400mM',...
    'LineStyle','none',...
    'Color',[1 0 0]);
set(errorbar1(1),'MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[0 0 0],...
    'DisplayName','0mM',...
        'LineStyle','none',...
    'Color',[0 0 0]);

% Create xlabel
xlabel('Tap','FontSize',24);

% Create ylabel
ylabel(msrname,'FontSize',24);

% Create legend
legend1 = legend(axes1,'show');
set(legend1,'EdgeColor',[1 1 1],'Location','EastOutside','YColor',[1 1 1],...
    'XColor',[1 1 1]);

