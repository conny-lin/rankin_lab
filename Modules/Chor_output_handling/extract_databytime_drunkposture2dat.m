function [DataTimed,errorMsg] = extract_databytime_drunkposture2dat(Data,time1,time2)
% downstream processing from import_shanespark_dat data

%% outputs
errorMsg = {};

%% legend
L = {...
'mwtnumber';
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

validateattributes(Data,{'table'},{'nonempty','ncols',datacolnumber},mfilename,'Data',1)

%% get tap data
i = Data.time >=time1 & Data.time <= time2;
DataTimed = Data(i,:);


end