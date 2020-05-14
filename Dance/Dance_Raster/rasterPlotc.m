function [figure1,Img, rastercolormap] = rasterPlotc(Data,time,visible,varargin)

%% DEFAULTS 
% graphics
axeslabel = false;
savelegend = false;
w = 3;
h = 3;
% raster color 
cmap = 'white';
grad = 'CL';
% define max min
gradmax = 0.6;
gradmin = -0.8;
% vararginprocessor
vararginProcessor
% ----------------------


%% process color map
switch cmap
    case 'white'
        % pause = white
        Blue = [0 0 1];
        Cyan = [0 1 1];
        White = [1 1 1];
        Yellow = [1 1 0];
        Red = [1 0 0];
        Dblue = [0 0 0.5];
        Dred = [0.5 0 0];
        step1 = 20;
        step2 = 60;
        Grad2 = [linspace(Blue(1),Cyan(1),step1+step2);linspace(Blue(2),Cyan(2),step1+step2);linspace(Blue(3),Cyan(3),step1+step2)];
        Grad3 = [linspace(Yellow(1),Red(1),step2);linspace(Yellow(2),Red(2),step2);linspace(Yellow(3),Red(3),step2)];
        Grad4 = [linspace(Red(1),Dred(1),step1);linspace(Red(2),Dred(2),step1);linspace(Red(3),Dred(3),step1)];
        rastercolormap =[Grad2'; White; Grad3'; Grad4'];
    case 'jet'
        rastercolormap = colormap('jet');
    case 'white intense'
        Blue = [0 0 1];
        Cyan = [0 1 1];
        White = [1 1 1];
        Yellow = [1 1 0];
        Red = [1 0 0];
        Dblue = [0 0 0.5];
        Dred = [0.5 0 0];
        step1 = 10;
        step2 = 30;
        Grad2 = [linspace(Blue(1),Cyan(1),step1+step2);linspace(Blue(2),Cyan(2),step1+step2);linspace(Blue(3),Cyan(3),step1+step2)];
        Grad3 = [linspace(Yellow(1),Red(1),step2);linspace(Yellow(2),Red(2),step2);linspace(Yellow(3),Red(3),step2)];
        Grad4 = [linspace(Red(1),Dred(1),step1);linspace(Red(2),Dred(2),step1);linspace(Red(3),Dred(3),step1)];
        rastercolormap =[Grad2'; White; Grad3'; Grad4']; 
        
end


%% define max color gradient
switch grad
    case 'CL'
        % check max min 3SD speed
        d = reshape(Data,numel(Data),1);
        f = d(d > 0);
%         speed_max = mean(f) + 3*std(f);
%         if speed_max > gradmax; warning('forward+3SD speed > %.1f',gradmax); end
        r = d(d < 0);
%         speed_min = mean(r) - 3*std(r);
%         if speed_min < gradmin; warning('reverse-3SD speed < %.1f',gradmin); end
        
        % convert speed 2 color based on gradmax
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

    case 'Evan'
        % forward movement
        For = Data;
        For(For<0)=0; % zeros for reversal
        For = For + 0.2; % add 0.2
        For(For == 0.2)= 0;

        % backwards
        Bac = Data;
        Bac(Bac>0)=0;
        Bac = Bac - 0.4;
        Bac(Bac==(-0.4))=0;

        % Tog
        Tog = For + Bac;
        Tog(Tog>0.8)=0.8; % max out at 0.8
        Tog(Tog<-0.8)=-0.8; % max out at -0.8
        Tog(Tog==0)=0.2; % if zero then 0.2
        Tog(1,1)=0.8; 
        Tog(1,2)=-0.8;
        
        Speedcolor = Tog;
end




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
set(figure1, 'PaperPosition',[0 0 w h],'PaperUnits','inches')

colormap(rastercolormap);
axes1 = axes('Parent',figure1);
hold(axes1,'on');

% set axes
% axes1 = axes('Parent',figure1,...
%     'YColor',[1 1 1],...
%     'XColor',[1 1 1],...
%     'Layer','top');
% %     'XTickLabel',timelabel,...
% %     'XTick',timeticks,...
% 

% expand graph to fill the page
sz = size(Speedcolor);
xlim([0 sz(2)])
ylim([0 sz(1)])


% remove axes
if ~axeslabel
    set(axes1,'Layer','top','XColor','none','YColor','none');
end

imagesc(Img,'Parent',axes1)


%% create raster legend------------------------------------
if savelegend
close all;
    figure1 = create_colorlegend(gradmin,gradmax,rastercolormap);

end
% -----------------------------------------------------------        




































