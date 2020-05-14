function [plateSumm,savename,MWTsourceID,recordTime,D] = rasterPlot_prepPlateSumm(path,start,finish,varargin)
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
frameInt = 0.2;
NWORMS = 200;
InputType = 2;
visibleG = 0;
saveoption = 0; % do not save trinity output
cleanup = 1;
%% VARARGIN PROCESSOR
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

%% remove
% %% survey folder content according to input type
% switch InputType
%     case 1 % group folder
%         [~,~,~,pMWT] = dircontent(path);
%     case 2 % mwt paths with the same group folder
%         pMWT = path;
%         [~,f] = cellfun(@fileparts,cellfun(@fileparts,pMWT,'UniformOutput',0),'UniformOutput',0);
%         if numel(unique(f)) ~= 1; 
%             error('mwt paths not all the same group'); 
%         end
% end

% %% IMPORT RAW DATA AND COMPILE TO MAT
% convertTrinityDat2Mat(pMWT);


%% declare matrix
expectedFrameN = (finish-start)/frameInt;
plateSumm =[];
recordTime = [];
MWTsourceID = [];
D = table;
D.mwtpath = pMWT;
D.data = cell(numel(pMWT),1);
for x = 1:numel(pMWT)
    fprintf('prcessing %d/%d MWT files\n',x,numel(pMWT));
    pMWTp = pMWT{x};
    [~,mwtn] = fileparts(pMWTp);
    [S,~,assayTime,d] = trinitySummary(pMWTp,start,finish,'frameInt',frameInt,'saveoption',saveoption);
    D.data{x} = d;
    if size(S,2) > expectedFrameN % if more than expected frameN
        warning('trinitySummary gave bigger frame N, delete extra [%s]',mwtn);
        S = S(:,1:expectedFrameN);
    elseif size(S,2) < expectedFrameN
        warning('trinitySummary gave smaller frame N, skip plate [%s]',mwtn);
    else
        plateSumm = [plateSumm;S];
        n = size(S,1);
        MWTsourceID = [MWTsourceID;repmat(x,n,1)];
        recordTime = [recordTime;repmat(assayTime,n,1)];
    end
end

% Evan's code-------------------------------------------------------------
% plateSumm = zeros((finish-start)*5+1,1)';
% groupResultSumm = [start:finish]'; 
% platecount = 0; %counting counts the group folders
% plateFrameTime = [start:frameInt:finish]'; % list frame time
% %Loops through folders
% for z = 1:numel(pMWT);
%     D = load([pMWT{z},'/trinitySummary.mat']);
%     D = D.masterData;
%     wormcount=1;
%     for wormi = 1:size(D,1);
%         % get data from this worm
%         storedData = D{wormi,2};
%         % if current trinity file contains data within start-finish
%         if storedData(1,1) < start && storedData(end,1) > finish;
%             wormcount=1+wormcount; % add count to k
%             % calculate data
%             storedData(storedData(:,1)<start | storedData(:,1)>finish,:)=[]; % delete data outside of start-finish
%             frameTimes = frameInt.* round(storedData(:,1)./.2); % round frame time to 0.2s increments
%             speed = storedData(:,4) .* storedData(:,6); % calculate speed by speed(s) in 4th column x bias (b) in 6th column
%             % average speed per frameTimes
%             frameTimesU = unique(frameTimes);
%             movemeans = zeros(size(frameTimesU));
%             movedir  = false(size(frameTimesU)); 
%             for jj = 1:numel(frameTimesU) 
%                 rowIndx = frameTimes==frameTimesU(jj);
%                 if abs(sum(speed(rowIndx,1))) == sum(abs(speed(rowIndx,1))); % if all movements within this frame range is in the same direction
%                     rowN = sum(rowIndx);
%                     movedir(jj)  = rowN>1; 
%                     movemeans(jj) = sum(speed(rowIndx,1))/rowN;
%                 else % if movement are different directions within frame
%                     conN = speed(rowIndx,1)<0; % get reversal 
%                     movemeans(jj) = sum(speed(conN,1))/sum(conN); % only calculate reversal mean
%                 end
%             end
%             % include data only if size of the means matches plate columns
%             if size(movemeans,1) == size(plateSumm,2);
%                 plateSumm = [plateSumm;movemeans'];
%             end
%         end
%     end
%     platecount = platecount+1; % increase plate count by one
% end
% plateSumm(any(isnan(plateSumm), 2), :) = []; % remove invalid data




%% save
% create save name
starttxt = num2str(start);
finishtxt =  num2str(round(finish));
ntext = num2str(size(Tog,1));
savename = ['rasterPlot_',starttxt,'_',finishtxt,'_N_',ntext];



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
