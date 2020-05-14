function OrgNewMWTintoGroupFolder(pM)
%% ORganize new MWT files
%% SEPARATE FILES
% pGF = '/Volumes/FLAME/Conny_Lin/20150910C_CL_100s30x10s10s_slo1rescue1';
% pM = '/Volumes/FLAME/Conny_Lin/20150910_100s30x10s10s_slo1rescue1';

% add path
addpath('/Users/connylin/Dropbox/MATLAB/Function_Library_Public');

%% GROUP CODES
% import
[~, ~, groupcode] = xlsread([pM,'/groupcode.xlsx'],'groupcode');
groupcode(cellfun(@(x) ~isempty(x) && isnumeric(x) && isnan(x),groupcode)) = {''};
% validate - no repeat codes or group name
if numel(unique(groupcode)) ~= numel(groupcode)
    error('group code or group name not unique')
else
    fprintf('PASS: group name and code are all unique\n')
end


%% PLATE CODE AND TIME
[~, ~, platetime] = xlsread([pM,'/groupcode.xlsx'],'MWTindex');
legend = platetime(1,:);
platetime = platetime(2:end,:);

% get tracker info
trackerlist = platetime(:,ismember(legend,'tracker'));
% get tracker name from exp name
[~,expname] =  fileparts(pM);
tracker_exp = regexp(expname,'(?<=^\d{8})[A-Z]{1}','match');
% find tracker not in current tracker
i = cellfun(@isempty,regexpi(trackerlist,tracker_exp));
platetime(i,:) = []; % remove info not in current tracker

% make third column = time
itime = find(ismember(legend,'start_time'));
% make first column = code
icode = find(ismember(legend,'platecode'));
platetime = platetime(:,[icode 2 itime]);


% platetime = platetime(2:end,1:3);
platetime(cellfun(@(x) ~isempty(x) && isnumeric(x) && isnan(x),platetime)) = {''};
platetime(cellfun(@isempty,platetime(:,1)),:)  = [];

% validate data - all MWT time must be unique
a = cellfun(@num2str,platetime(:,3),'UniformOutput',0);
errormsg = 'plate time not unique';
passmsg = 'PASS: plate times are all unique\n';
if numel(unique(a)) ~= numel(a)
    error(errormsg)
else
    fprintf(passmsg)
end
% validate data - all MWT time must be unique
a = platetime(:,1);
errormsg = 'plate code not unique';
passmsg = 'PASS: plate code are all unique\n';
if numel(unique(a)) ~= numel(a)
    % find repeat names
    b = tabulate(a);
    i = cell2mat(b(:,2)) > 1;
    display 'repeated codes:';
    display(platetime(i,1))  
    error(errormsg)
else
    fprintf(passmsg)
end


%% sort by plate time
platetime = sortrows(platetime,find(strcmp(legend,'start_time')));


%% match plate number with plate time entries
% get MWT plate time stamp
[~,~,MF,pMFlist] = dircontent(pM);
a = regexpcellout(MF,'_','split');
timestamp = zeros(size(a,1),1);
for x = 1:size(a,1)
    timestamp(x,1) = str2num(a{x,2}(1,1:4));
end

%% match time
PT = platetime(:,[1 3]);
for x = 1:size(PT,1);
    t = PT{x,2};
    i = ismember(timestamp,t);
    if sum(i) == 1
        PT{x,3} = timestamp(i);
        PT(x,4) = MF(i);
        PT(x,5) = pMFlist(i);
        timestamp(i) = 0; % remove taken ones

    end
end

%% see missing ones
i = cellfun(@isempty,PT(:,3));
j = find(i);
for x = 1:numel(j)
    
    % display before and after
    % if the first one
    if j(x) == 1
        disp(PT(j(x):j(x)+1,1:3))
    elseif j(x) == size(PT,1)
        disp(PT(j(x)-1:j(x),1:3))
    else
        disp(PT(j(x)-1:j(x)+1,1:3))
    end
    opt = cellstr(num2str((timestamp(timestamp ~= 0))));
    a = chooseoption([opt;{'none of the above'}],1,'Enter 0 is nothing matches');
    if strcmp(a,'none of the above') ~= 1
       
        i = timestamp == str2num(a);
        
        PT{j(x),3} = timestamp(i);
        PT(j(x),4) = MF(i);
        PT(j(x),5) = pMFlist(i);
        timestamp(i) = 0; % remove taken ones
    end
end

%% check MWT without assignment
disp(PT(:,1:3));
display 'MWT without assignment:'
disp(MF(timestamp ~= 0))


%% ok?
if input('is that ok? (y=1 n=0) ') == 0
    disp('manually change fix PT matrix then run folder change')
end

%% move folder
% remove ones without path
i = cellfun(@isempty,PT(:,3));
PTM = PT(~i,:);
display('starting to move');
% get folder name
for x = 1:size(PTM,1)
   foldername = groupcode{ismember(groupcode(:,1),PT{x,1}(1)),2};
   PT(x,6) = {foldername};
   if isempty(PT(x,5)) == 0
       fT = [pM,'/',foldername];
       if isdir(fT) == 0; mkdir(fT); end
       s = PTM{x,5};
       d = fileparts(regexprep(s,fileparts(s),fT));
       movefile(s,d,'f');
   end
end


%% create report

% platetime(:,1:2)
T = cell2table(PT(:,[1 6 4]),'VariableNames',{'platecode','groupname','MWTname'});
cd(pM);
writetable(T,'plateinfo.csv');

%% 
display 'FINISHED';






















