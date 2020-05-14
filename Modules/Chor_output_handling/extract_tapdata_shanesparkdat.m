function [TapData,errorMsg] = extract_tapdata_shanesparkdat(Data)
% downstream processing from import_shanespark_dat data

%% outputs
errorMsg = {};

% legend
L = {...
'mwtnumber';
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

validateattributes(Data,{'table'},{'nonempty','ncols',datacolnumber},mfilename,'Data',1)

% get tap data
i = Data.tap == 1;
TapData = Data(i,:);


end