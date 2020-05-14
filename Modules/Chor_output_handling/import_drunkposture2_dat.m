function [Data,errorMsg,LM] = import_drunkposture2_dat(pMWT,outputtype)
%% chor script
%% odrunkposture2 = '-O drunkposture2 -o nNslwakbcemM';
% pMWT can be paths to actual file or path to mwt folder

%% default
if nargin < 2
   outputtype = 'table'; 
end

%% CHECK INPUT VARIABLES
i = cellfun(@isdir,pMWT);
if sum(i) > 0
    fileinputtype = 'folder'; 
else
    fileinputtype = 'datfile';
end

%% legend
L = {...
'time'; % -- always the first column unless included again
'ntracked'; % -- the number of objects tracked 
'goodnumber'; %  -- the number of objects passing the criteria given
'speed'; %  -- speed of movement
'length';% -- measured along major axis, not curve of object ';
'width';% -- width of the rectangle framing the body';
'aspect';% length/width ';
'kink';% -- head/tail angle difference from body (in degrees) ';
'bias'; % fractional excess of time spent moving one way ';
'curve'; % -- average angle (in degrees) between body split into 5 segments ';
'area'; % body area';
'midline';% -- length measured along the curve of object ';
'morphwidth' % -- mean width of body about midline ';
}; 

datacolnumber = numel(L);
extname = 'drunkposture2.dat';
LM = [{'mwtid'};L];


%% set up 
Data = cell(size(pMWT,1),1);
npmwt = numel(pMWT);
errorMsg = cell(size(pMWT));


%% get files
for m = 1:npmwt    
    % get path
    if strcmp(fileinputtype,'folder')
        pmwt = pMWT{m}; 
        [~,p] = dircontent(pmwt,['*.',extname]); 
        % delete temperarary .*. files
        if size(p,1) ~= 1
            cellfun(@delete,p(regexpcellout(p,['\<[.]\w{1,}[.]',extname,'\>'])))
            [~,p] = dircontent(pmwt,['*.',extname]);
        end
        % if no file found
        if isempty(p) == 1; errorMsg{m} = 'no file'; continue; end
        % if more than one file found
        if numel(p) > 1
           p = p(1); % use the first file
           errorMsg{m} = 'multiple files';
        end
        pmwt = char(p);
    else
        pmwt = pMWT{m};  
    end
    % import
    d = dlmread(pmwt);
    % check import file size
    if size(d,2) ~= datacolnumber
        errorMsg{m} = 'file col number incorrect'; continue
    else
        % create mwt number
        mn = repmat(m,size(d,1),1);
        d = [mn d];
        if strcmp(outputtype,'table')
            d = array2table(d,'VariableNames',LM);
        end
        Data{m} = d;
    end
end     

    


