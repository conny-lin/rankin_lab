function figure1 = rasterPlot_colorSpeed_hot(plateSumm,varargin)

%% default
NWORMS = Inf; %NWORMS = 200;
visibleG = 0;
vararginProcessor;



%% summarizing
% forward movement
For = plateSumm;
For(For<0)=0; % zeros for reversal
For = For + 0.2; % add 0.2
For(For == 0.2)= 0;

% backwards
Bac = plateSumm;
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
 

%% scale output to N = scaleNumWorms
nWorms = size(Tog,1);
if isinf(NWORMS) == 0
    if nWorms < NWORMS % if less worm than defined in scaleNumWorms, add zeros
        extr = NWORMS - nWorms;
        extrRows = zeros(extr, size(Tog,2));
        Tog = [extrRows;Tog]; 
    else % take the first 200 worms
        Tog(NWORMS:end,:)=[];
    end
end



%% create image
if visibleG == 0
    figure1 = figure('Color',[1 1 1],'Visible','off');
else
    figure1 = figure('Color',[1 1 1],'Visible','on');
end
colormap jet;
axes1 = axes('Parent',figure1,...
    'ZColor',[0.8 0.8 0.8],...
    'YDir','reverse',...
    'YColor',[0.8 0.8 0.8],...
    'XColor',[0.8 0.8 0.8],...
    'Layer','top');
imagesc(Tog,'Parent',axes1)