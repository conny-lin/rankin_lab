function D = transform_roundSpeedbytime(D,varargin)
%% transform_rounddatabytime(D)
% D must be a table and has time
% last updated: 201603061028

%% varargin
overwriteTime = 0;
timefactor = .1;
%% process varargin
vararginProcessor

%%
D.timeround = round(D.time.*(timefactor*100))./(timefactor*100);
a = tabulate(D.timeround);
a(a(:,2)<2,:) = [];
tru = a(:,1);
for trui = 1:numel(tru)
    s = []; b = []; % reset critical values
    imulti = find(D.timeround==tru(trui)); % get time within time frame
    bias = D.bias(imulti); % get bias
    speed = D.speed(imulti); % get speed
    bu = unique(D.bias(imulti)); % get unique bias
    % deal with bias containing nan value
    if any(isnan(bu)) % remove NaN bias
        bu(isnan(bu)) = [];
    end
    if isempty(bu) % if no bias left after removing nan, 
        D.speed(imulti(1)) = nanmean(speed); % calculate speed without nan data
        D.bias(imulti(1)) = NaN; % store bias as nan
        D(imulti(2:end),:) = []; % remove duplicated data point
        continue
    end
    % deal with bias do not contain nan value
    if numel(bu)==1 % if only one bias
        s = nanmean(D.speed(imulti)); % take mean
        b = bu; % mark bias
    elseif numel(bu)>1 % if more than one bias
        if ismember(bu,[0 -1]) % if pause and reverse, prefer reverse
           s = nanmean(speed(bias~=0));
           b = -1;
        elseif ismember(bu,[0 1]) % if pause and forward, prefer fowrward
            s = nanmean(speed(bias~=0));
            b = 1;
        elseif ismember(bu,[-1 1]) % if pause or forward, 
            st = nanmean(speed.*bias); % take a mean of s
            if st>0 % if mean is more positive
                s = nanmean(speed(bias==1)); % take only mean from fowarward
                b=1; % mark it as forward speed
            elseif st<0 % if mean is more negative
                s = nanmean(speed(bias==-1)); % take only mean from reverse
                b = -1; % mark it as reverse speed
            elseif st==0 % if mean equal to zero
                s = 0; % mark speed as zero
                b=0; % mark bias as zero
            else
                error('code for this scenario');
            end
        elseif ismember(bu,[0 -1 1]) % if all 3 move are included
            s = nanmean(speed.*bias); % take mean of all movement
            b = 0; % mark it as pause
        else
            error('stop to process');
        end
    else
        error('code for this scenario');
    end
    % update data
    D.bias(imulti(1)) = b; % record bias
    D.speed(imulti(1)) = s; % record speed
    D(imulti(2:end),:) = []; % delete the rest
end

%% overwrite time
if overwriteTime
    % store new time and averaged data
    D.time = D.timeround;
    D.timeround= [];
end
            
            
%% validate
if numel(unique(D.time)) ~= numel(D.time)
    error('failed');
end