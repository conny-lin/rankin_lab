function [figure1,Tog] = rasterPlot_color(plateSumm,varargin)
%% RL201510211310_rastor_plot create rastor plot from a MWT analyzed folder 
% [pSave,Motherfoldercontents] = RL201510211310_rastor_plot(pMWT,start,finish,pSave)
% path is path to an folder with muliple MWT plate folder in string start is 
% start time of analysis in seconds finish is end time of analysis in seconds
% the frames will be aligned to tap, if tap exist
%
%     Examples:
%         start = 80;
%         finish = 100;
%     options:
%         'frameInt',0.2 = setting frame interval to every 0.2seconds
%         'NWORMS' = a numeric value setting max number of worm 
%             included in the output, or inf to set output to all worms
%         'InputType': 
%             1 = path to one group folder (default), delete trinity.*.dat
%             after creating trinitySummary.mat
%             2 = path to one multiple MWT of the same group; if not the same
%             group will return an error; 
%         saveoption = 0, do not save trinity summary in database
%         cleanup = 1, delete trinity.dat file after summarizing it to .mat
%
%     OUTPUT:
%         output frame time is rounded to 0.2s increments
% 
%     REQUIREMENT 
%         .evanall.dat (oevanall = '-O evanall -N all -o nNss*b12')
%         .trinity.dat -o tnNss*b12xyMmeSakcr'; 
% 
%     MODIFIED from speed_analysis_code from Alex 20151021



%% DEFAULTS
nInput = 3;
% frameInt = 0.2;
NWORMS = Inf;
% InputType = 2;
visibleG = 0;
% saveoption = 0; % do not save trinity output
% cleanup = 1;
%%
vararginProcessor


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
if ~visibleG
    figure1 = figure('Color',[1 1 1],'Visible','off');
else
    figure1 = figure('Color',[1 1 1],'Visible','on');
end
colormap('jet');
axes1 = axes('Parent',figure1,...
    'ZColor',[0.8 0.8 0.8],...
    'YDir','reverse',...
    'YColor',[0.8 0.8 0.8],...
    'XColor',[0.8 0.8 0.8],...
    'Layer','top');
imagesc(Tog,'Parent',axes1)


%% remove
% save
% % create save name
% starttxt = num2str(start);
% finishtxt =  num2str(round(finish));
% ntext = num2str(size(Tog,1));
% savename = ['rasterPlot_',starttxt,'_',finishtxt,'_N_',ntext];



end









%% SUBFUNCTIONS
function [varargout] = dircontent(p,varargin)
% a = full dir, can be folder or files, b = path to all files, 
% c = only folders, d = paths to only folders

%% get dir
switch nargin
    case 0
        varargout{1} = [];
        return
    
    case 1
        % get dir
        [a,b,c,d] = dirc(p);
        varargout{1} = a;
        varargout{2} = b;
        varargout{3} = c;
        varargout{4} = d;
        return
    
    
    case 2
        ext = varargin{1};
        % dircontnetext
        cd(p); % go to directory
        a = {};
        a = dir(ext); % list content
        a = {a.name}'; % extract folder names only
        a(ismember(a,{'.','..','.DS_Store'})) = []; 
        b = {};
        for x = 1:size(a,1); % for all files 
            b{x,1} = [p,'/',a{x,1}]; % make path for files
        end
         varargout{1} = a;
         varargout{2} = b;
         return
    
         
    case 3
        if strcmp(varargin{1},'Option')
            option = varargin{2};
        end
        
end


switch option
    case 'MWT'
        p1 = regexp(genpath(p),':','split')';
        [~,fn2] = cellfun(@fileparts,p1,'UniformOutput',0);
        i = ismember(fn2,{''});
        fn2(i) = []; 
        p1(i) = [];
        i = regexpcellout(fn2,'\<\d{8}[_]\d{6}\>');
        varargout{1} = fn2(i);
        varargout{2} = p1(i);
    case 'MWTall'
        p1 = regexp(genpath(p),':','split')';
        [~,fn2] = cellfun(@fileparts,p1,'UniformOutput',0);
        i = ismember(fn2,{''});
        fn2(i) = []; 
        p1(i) = [];
        i = regexpcellout(fn2,'\<\d{8}[_]\d{6}');
        varargout{1} = fn2(i);
        varargout{2} = p1(i);

    otherwise
        error 'No such option';
end
end




%% Subfunction
function [a,b,c,d] = dirc(p)
% a = full dir, can be folder or files, b = path to all files, 
% c = only folders, d = paths to only folders
cd(p); % go to directory
a = {}; % create cell array for output
a = dir; % list content
a = {a.name}'; % extract folder names only
a(ismember(a,{'.','..','.DS_Store'})) = []; 
b = {};
c = {};
d = {};
for x = 1:size(a,1); % for all files 
    b(x,1) = {strcat(p,'/',a{x,1})}; % make path for files
    if isdir(b{x,1}) ==1; % if a path is a folder
        c(end+1,1) = a(x,1); % record name in cell array b
        d(end+1,1) = b(x,1); % create paths for folders
    else
    end
end
end
