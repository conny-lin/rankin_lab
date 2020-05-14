function [DataOut,legend] = extractChorData(Data,colName,varargin)
% downstream processing from import_shanespark_dat data

%% addpath
pFun = {'/Users/connylin/Dropbox/Code/Matlab/Library/General'};
for x = 1:numel(pFun); addpath(pFun{x}); end
 
%% default
time1 = 0;
time2 = Inf;
outputtype = 'table';
% % legend
% colName = {...
% 'mwtnumber';
% 'time'; % -- always the first column unless included again
% 'ntracked'; % -- the number of objects tracked 
% 'goodnumber'; %  -- the number of objects passing the criteria given
% 'speed'; %  -- speed of movement
% 'length';% -- measured along major axis, not curve of object ';
% 'width';% -- width of the rectangle framing the body';
% 'aspect';% length/width ';
% 'kink';% -- head/tail angle difference from body (in degrees) ';
% 'bias'; % fractional excess of time spent moving one way ';
% 'curve'; % -- average angle (in degrees) between body split into 5 segments ';
% 'area'; % body area';
% 'midline';% -- length measured along the curve of object ';
% 'morphwidth' % -- mean width of body about midline ';
% }; 
varnames = colName;
vararginProcessor;

%% process input
datacolnumber = numel(colName);
if istable(Data) == 0
    Data = array2table(Data,'VariableNames',colName);
end
validateattributes(Data,{'table'},{'nonempty','ncols',datacolnumber},mfilename,'Data',1)
% get time data
i = Data.time >=time1 & Data.time <= time2;
D = Data(i,:);
% get specific col
v = D.Properties.VariableNames;
i = ismember(v,varnames);
D(:,~i) = [];
legend = D.Properties.VariableNames;
% process output type
switch outputtype
    case 'table'
        DataOut = D;
    case 'array'
        DataOut = table2array(D);
end

end