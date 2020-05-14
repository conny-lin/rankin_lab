%%%%%%for ALL Plate
function plateSumm = RL201510211310_rastor_plot(pG,start,finish,varargin)
%% RL201510211310_rastor_plot create rastor plot from a MWT analyzed folder 
% [pSave,Motherfoldercontents] = RL201510211310_rastor_plot(pMWT,start,finish,pSave)
%     pG is path to an folder with muliple MWT plate folder in string
%     start is start time of analysis in seconds
%     finish is end time of analysis in seconds
%     Examples:
%         start = 80;
%         finish = 100;
%     options:
%         'frameInt',0.2 = setting frame interval to every 0.2seconds
%         'NWORMS' = a numeric value setting max number of worm 
%             included in the output, or inf to set output to all worms
%
%     output:
%     output frame time is rounded to 0.2s increments
% 
%     REQUIREMENT 
%     .evanall.dat (oevanall = '-O evanall -N all -o nNss*b12')
%     .trinity.dat -o tnNss*b12xyMmeSakcr'; 
% 
%     MODIFIED from speed_analysis_code from Alex 20151021



%% default settings
frameInt = 0.2;
NWORMS = 200;

if nargin > 3
    A = varargin;
    
    if isinteger(numel(A)/2) == 1
        error('variable entries incorrect')
    end
    callName = A([1:2:numel(A)]);
    setValue = cell2mat(A([2:2:numel(A)]));
    for x = 1:numel(callName)
        str = sprintf('%s = %d;',callName{x},setValue(x));
        eval(str)
    end

end

%% survey folder content
[~,~,~,pMWT] = dircontent(pG);

%% transform data in mat file
for z = 1:numel(pMWT); 
    [~,TrinityD] = dircontent(pMWT{z},'trinitySummary.mat');
    if numel(TrinityD) == 0
        [~,pTrinity] = dircontent(pMWT{z},'*.trinity.*.dat');
        masterData = cell(size(pTrinity,1),2); % declare output array
        % get worm ID
        [~,f] = cellfun(@fileparts,pTrinity,'UniformOutput',0);
        % record worm id
        masterData(:,1) = regexpcellout(f,'(?<=trinity[.])\d{1,}','match');
        for wormi = 1:numel(pTrinity);
            masterData{wormi,2} = dlmread(pTrinity{wormi}); % read .trinity.*.dat data from a worm
        end
        save([pMWT{z},'/trinitySummary.mat'],'masterData','-v7.3');
    end
end


%% declare matrix
plateSumm = zeros((finish-start)*5+1,1)';
groupResultSumm = [start:finish]';
platecount = 0; %counting counts the group folders
% list frame time
plateFrameTime = [start:frameInt:finish]';

%Loops through folders
for z = 1:numel(pMWT);
    D = load([pMWT{z},'/trinitySummary.mat']);
    D = D.masterData;
%     [~,TrinityD] = dircontent(pMWT{2},'trinitySummary.mat');
%     [~,pTrinity] = dircontent(pMWT{z},'*.trinity.*.dat');
    
    wormcount=1;
    for wormi = 1:size(D,1);
%         storedData = dlmread(TrinityD{wormi}); % read .trinity.*.dat data from a worm
        storedData = D{wormi,2};
        % if current trinity file contains data within start-finish
        if storedData(1,1) < start && storedData(end,1) > finish;
            wormcount=1+wormcount; % add count to k

            storedData(storedData(:,1)<start | storedData(:,1)>finish,:)=[]; % delete data outside of start-finish
            frameTimes = frameInt.* round(storedData(:,1)./.2); % round frame time to 0.2s increments
            speed = storedData(:,4) .* storedData(:,6); % calculate speed by speed(s) in 4th column x bias (b) in 6th column

            % average speed per frameTimes
            frameTimesU = unique(frameTimes);
            movemeans = zeros(size(frameTimesU));
            movedir  = false(size(frameTimesU)); 
            for jj = 1:numel(frameTimesU) 

                rowIndx = frameTimes==frameTimesU(jj);
                if abs(sum(speed(rowIndx,1))) == sum(abs(speed(rowIndx,1))); % if all movements within this frame range is in the same direction
                    rowN = sum(rowIndx);
                    movedir(jj)  = rowN>1; 
                    movemeans(jj) = sum(speed(rowIndx,1))/rowN;
                else % if movement are different directions within frame
                    conN = speed(rowIndx,1)<0; % get reversal 
                    movemeans(jj) = sum(speed(conN,1))/sum(conN); % only calculate reversal mean
                end

            end
            % include data only if size of the means matches plate columns
            if size(movemeans,1) == size(plateSumm,2);
                plateSumm = [plateSumm;movemeans'];
            end
        end
    end

    platecount = platecount+1; % increase plate count by one

end


%% summarizing
plateSumm(any(isnan(plateSumm), 2), :) = []; % remove invalid data

% forward movement
For = plateSumm;
For(For<0)=0; % zeros for reversal
For = For + 0.2; % add 0.2
For(For == 0.2)=0;

% backwards
Bac = plateSumm;
Bac(Bac>0)=0;
Bac = Bac - 0.4;
Bac(Bac==(-0.4))=0;

% Tog
Tog = For + Bac;
Tog(Tog>0.8)=0.8;
Tog(Tog<-0.8)=-0.8;
Tog(Tog==0)=0.2;
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
imagesc(Tog)


%% save data
starttxt = num2str(start);
finishtxt =  num2str(round(finish));
ntext = num2str(size(Tog,1));
savefigeps(['rasterPlot_',starttxt,'_',finishtxt,'_N_',ntext],pG); % save data to group folder
save([pG,'/',starttxt,'_',finishtxt,'_N_',ntext,'.rasterData'],'plateSumm','-ascii');
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
