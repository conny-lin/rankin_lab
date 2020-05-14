function [figure1,axes1] = Graph_errorbar(X,Y,E,gname,colorset,varargin)
%% Graph_errorbar
% based on makefig_Errorbar_line

%% DEFAULTS & VARARGIN
nInput = 5;
visiblesetting =  'on';
xname = 'Tap';
yname = '';
titlestr = '';
% varargin processer
if nargin > nInput
    A = varargin;
    if isinteger(numel(A)/2) == 1
        error('variable entries incorrect')
    end
    callName = A(1:2:numel(A));
    % check if call names are strings
    for x = 1:numel(callName)
       if ischar(callName{x}) == 0
           error('varargin incorrect');
       end
    end
    % get values
    setValue = A(2:2:numel(A));
    for x = 1:numel(callName)
        % determine input type
        if eval(sprintf('ischar(%s)',callName{x})) == 1
            eval(sprintf('%s = ''%s'';',callName{x},setValue{x}));
        elseif eval(sprintf('isnumeric(%s)',callName{x}))
            eval(sprintf('%s = %d;',callName{x},cell2mat(setValue(x))));
        elseif eval(sprintf('iscell(%s)',callName{x}))
            a = setValue{x};
            eval(sprintf('%s = a;',callName{x}));
        else
            error('need coding');
        end
    end
end

% variables post-processing
xname = regexprep(xname,'_',' ');
yname = regexprep(yname,'_',' ');
gname = regexprep(gname,'_',' ');


%% Create figure
figure1 = figure('Color',[1 1 1],'Visible',visiblesetting);

axes1 = axes('Parent',figure1,'FontSize',18);
box(axes1,'off');

xmin = min(min(X))*0.9;
xmax = max(max(X)*1.1);
if isnan(xmin) == 0 && isnan(xmax) == 0
    xlim(axes1,[xmin xmax]);
elseif isnan(xmin) == 1
    xlim(axes1,[0 xmax]);
elseif isnan(xmax) == 0;
end

% ymin = min(min(Y-E))*0.9;
% ymax = max(max(Y+E)*1.1);
% if isnan(ymin) == 0 && isnan(ymax) == 0
%     ylim(axes1,[ymin ymax]);
% elseif isnan(ymin) == 1
%     ylim(axes1,[0 ymax]);
% elseif isnan(ymax) == 0;
% end

hold(axes1,'all');
errorbar1 = errorbar(X,Y,E,'LineWidth',1);

n = min([size(colorset,1) size(Y,2)]);
for x = 1:n
    set(errorbar1(x),'Color',colorset(x,1:3))
    
end

for x = 1:size(Y,2)
    set(errorbar1(x),'DisplayName',gname{x})
end
title(titlestr,'FontSize',12);

xlabel(xname,'FontSize',18);

% ylable 
ylabel(yname,'FontSize',18);

% legend
legend1 = legend(axes1,'show');
set(legend1,'EdgeColor',[1 1 1],...
    'Location','NorthEastOutside','FontSize',12);

%% text strings (suspended r20151203)

% a = PlateN;
% b = '';
% for x = 1:numel(a)
%     b = [b,num2str(a(x)),', '];
% end
% 
% str1 = ['N(plate) = ',b(:,1:end-2)];
% 
% textboxstr = [str1];
% annotation(figure1,'textbox',...
%     [0.70 0.015 0.256 0.05],...
%     'String',textboxstr,...
%     'FitBoxToText','off',...
%     'EdgeColor',[1 1 1]);


