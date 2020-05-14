function [speed,bias,movedur,pretapbias,pretapspeed,D] = calculate_biasedSpeed_poststim(D,varargin)


% assaytype = 'poststim'; 
%     need to include first row as pre-tapped response, second row as tapped 
%     response, third+rows as responses for analysis
% assaytype pretim, 
%     include responses before the tap (exclude tap)

%% process varargin
responseNum = 2; % need to include first row as pre-tapped response, second row as tapped response, third+rows as responses for analysis
% assaytype pretim, include responses before the tap (exclude tap)
vararginProcessor;


%% calculation
movedur = repmat(Inf,responseNum,1);
speed = movedur; 
bias = movedur; 
pretapbias = Inf; pretapspeed = Inf;
%% get bias tap info
itap = find(D.tap==1);
pretapbias = D.bias(itap-1);
% if isnan(pretapbias); pretapbias=4; end
tapbias = D.bias(itap);
% if isnan(tapbias); tapbias=4; end
pretapspeed = D.speed(itap-1);
tapspeed = D.speed(itap);
if pretapbias~=tapbias; 
   PreTap = D(1:itap-1,:);
   D(1:itap-1,:) = []; 
else
   PreTap = D(1:itap,:);
   D(1:itap,:) = []; 
end
if isempty(D)==1; return; end
% find movement variations
df = find([false;diff(D.bias)]);
if numel(df)>0; 
    b = [D.bias(1);D.bias(df)]';
else
    b = D.bias(1);
end
% find movement directions
dir = unique(b); % unique directions
if any(isnan(dir))
    dir(isnan(dir)) = [];
    dir = [dir NaN];
end
% first movement bias
if isnan(D.bias(1));
    firstdir = NaN;
else
    firstdir = D.bias(1);
end

%% decide bias and valid speed data
if numel(dir)==1 % if all going one direction
    speed = nanmean(D.speed);
    bias=firstdir;
    movedur = D.time(end)-D.time(1);
elseif numel(dir)>1
    if df(end)~=size(D,1)
        dff = [1;df;size(D,1)];
    else
        dff=[1;df];
    end
    bias = b';
    speed = nan(size(bias));
    movedur = speed;
    for dfi = 1:numel(dff)-1
        s = nanmean(D.speed(dff(dfi):dff(dfi)+1));
        if b(dfi)==-1
            speed(dfi) =-s;
        elseif b(dfi)==0 || b(dfi)==1;
            speed(dfi) =s;
        elseif isnan(b(dfi))
            speed(dfi) =NaN;
        end
    end
    df2 = [1;df-1;size(D,1)];
    for dfi = 1:numel(df2)-1
        movedur(dfi) = D.time(df2(dfi+1))-D.time(df2(dfi));
    end
else
    error('code for this scenario');
end

%% deal with number of outputs
n = numel(speed);
if responseNum ==1
    a = [speed bias movedur];
    a(2:end,:) = [];
    speed = a(1); bias = a(2); movedur = a(3);
elseif responseNum > 1 
    if n < responseNum
        a = nan(responseNum,3);
        a(1:n,1) = speed;
        a(1:n,2) = bias;
        a(1:n,3) = movedur;
        speed = a(:,1); 
        bias = a(:,2);
        movedur = a(:,3);
    elseif n > responseNum
        speed(responseNum+1:end) = [];
        bias(responseNum+1:end) = [];
        movedur(responseNum+1:end) = [];
    end
end





















