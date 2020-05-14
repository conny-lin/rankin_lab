function [Data,errorMsg,LM] = import_shanespark_dat(pMWT)
%% chor script
% oshanespark = '-O shanespark -o nNss*b12M'; % standard for Beethoven


%% legend
L = {...
'time'; % -- always the first column unless included again
'ntracked'; % -- the number of objects tracked 
'goodnumber'; %  -- the number of objects passing the criteria given
'speed'; %  -- speed of movement
'speed_std'; %  -- speed of movement
'bias'; %  -- fractional excess of time spent moving one way 
'tap'; %  -- whether a tap (stimulus 1) has occurred 
'puff'; %  -- whether a puff (stimulus 2) has occurred
'morphwidth'}; % -- mean width of body about midline 
datacolnumber = numel(L);
extname = 'shanespark.dat';
LM = [{'mwtnumber'};L];

%% conversion
Data = cell(size(pMWT,1),1);
npmwt = numel(pMWT);
errorMsg = cell(size(pMWT));
for m = 1:npmwt    
    pmwt = pMWT{m};    
    % import files
    [~,mwtname] = fileparts(pmwt);
    mwtnamenumber = str2double(regexprep(mwtname,'_',''));    
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
    % import
    d = dlmread(char(p));
    % check import file size
    if size(d,2) ~= datacolnumber
        errorMsg{m} = 'file col number incorrect'; continue
    else
        % create mwt number
        mn = repmat(mwtnamenumber,size(d,1),1);
        d = [mn d];
        
        data = array2table(d,'VariableNames',LM);
        Data{m} = data;
    end
end


