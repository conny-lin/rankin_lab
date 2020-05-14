function D = transform_roundSpeedbytime2(D,varargin)
%% transform_rounddatabytime(D)
% D must be a table and has time
% last updated: 201603061028

%% varargin
overwriteTime = 0;
timefactor = .1;
%% process varargin
vararginProcessor

%% set paramters
D.timeround = round(D.time.*(timefactor*100))./(timefactor*100);
a = tabulate(D.timeround);
a(a(:,2)<2,:) = [];
tru = a(:,1);
if isempty(tru); return; end
%% construct evaluation matrix
TR = cell(numel(tru),5);
for trui = 1:numel(tru)
    s = []; b = []; % reset critical values
    TR{trui,1} = tru(trui);
    imulti = find(D.timeround==tru(trui)); % get time within time frame
    TR{trui,2} = imulti;
    TR{trui,3} = imulti(1); % output row
    bias = D.bias(imulti); % get bias
    speed = D.speed(imulti); % get speed
    bu = unique(D.bias(imulti)); % get unique bias
    if any(isnan(bu))
       bu(isnan(bu))=[]; 
       if ~isempty(bu); % if there are other data other than NaN bias
           % delete NaN data
           speed(isnan(bias)) = [];
           bias(isnan(bias)) = [];
       else
           bu=NaN;
       end
    end
    TR{trui,4} = bias;
    TR{trui,5} = speed;
    TR{trui,6} = bu;
end
TR = cell2table(TR,'VariableNames',{'time','rowid','rowid_out','bias','speed','biasu'});
% create output vectors
TR.speedout= Inf(size(TR,1),1);
TR.biasout= Inf(size(TR,1),1);
%% deal with bias containing nan value
if iscell(TR.biasu)==1
    nanVal = cellfun(@sum,cellfun(@isnan,TR.biasu,'UniformOutput',0))>0;
else
    nanVal = isnan(TR.biasu);
end
if sum(nanVal) > 0
    C = TR(nanVal,:);
    if iscell(C.biasu)==1
        i = cellfun(@isnan,C.biasu); % validate only one nan item
    else
        i = isnan(C.biasu); % validate only one nan item
    end
    TR.speedout(nanVal) = cellfun(@nanmean,C.speed);
    TR.biasout(nanVal) = nan(size(C,1),1);
end
%% get unique number count
if iscell(TR.biasu)==1
    buN = cellfun(@numel,TR.biasu);
else
    buN = ones(size(TR.biasu));
end
%% process bias == 1
if any(buN)==1
    C = TR(buN==1,:);
    if iscell(C.biasu);
        TR.biasout(buN==1) = cell2mat(C.biasu);
    else
        TR.biasout(buN==1) = C.biasu;
    end
    TR.speedout(buN==1) = cellfun(@nanmean,C.speed);
end
%% deal with bias unique > 1
if any(buN>1)
    %% deal with biasu with 2 directions
    %% for ones containig zero bias
    ind = find(buN==2);
    ind = ind(any(cell2mat(TR.biasu(ind)')==0));    
    if ~isempty(ind)
        C = TR(ind,:);
        r = cellfun(@numel,C.bias);
        if numel(unique(r))>1
            b = nan(max(r),size(C,1));
            s = b;
            for ci = 1:size(b,2)
               s(1:r(ci),ci) = C.speed{ci};
               b(1:r(ci),ci) = C.bias{ci};
            end
            b(b==0) = NaN;
        else
            b = cell2mat(C.bias');
            s = cell2mat(C.speed');
        end
        s(b==0) = NaN;
        TR.speedout(ind) = nanmean(s);
        b(b==0) = NaN;
        TR.biasout(ind) = nanmean(b);
    end
    %% for ones containig [-1 1]
    ind = find(buN==2);
    ind = ind(~any(cell2mat(TR.biasu(ind)')==0));
    if ~isempty(ind)
        C = TR(ind,:);
        %% get index
        i = [0;cumsum(cellfun(@numel,C.bias))];
        g = nan(max(i),1);
        for ii = 1:numel(i)-1
            g(i(ii)+1:i(ii+1))=ii;
        end
        %% calculate mean
        s = cell2mat(C.speed);
        b = cell2mat(C.bias);
        a = [g s b];
        [gn,sm] = grpstats(a(:,2).*a(:,3),a(:,1),{'gname','mean'});
        gn = cellfun(@str2num,gn);
        %% for positive mean
        b = a(ismember(a(:,1),gn(sm>0)),:);
        if ~isempty(b)
            b(b(:,3)<0,:) = [];
            [rowi,s2] = grpstats(b(:,2),b(:,1),{'gname','mean'});
            rowi = cellfun(@str2num,rowi);
            C.speedout(rowi) = s2;
            C.biasout(rowi) = 1;
        end
        %% for negative mean
        b = a(ismember(a(:,1),gn(sm<0)),:);
        if ~isempty(b)
            b(b(:,3)>0,:) = [];
            [rowi,s2] = grpstats(b(:,2),b(:,1),{'gname','mean'});
            rowi = cellfun(@str2num,rowi);
            C.speedout(rowi) = s2;
            C.biasout(rowi) = -1;
        end
        %% for zero mean
        b = a(ismember(a(:,1),gn(sm==0)),:);
        if ~isempty(b)
            C.speedout(gn(sm==0)) = 0;
            C.biasout(gn(sm==0)) = 0;
        end
        TR.speedout(ind)=C.speedout;
        TR.biasout(ind)=C.biasout;
    end   
    %% check for 2 or more types, average and assign sign
    
    ind = find(buN>2);
    if ~isempty(ind)
       C = TR(ind,:);
       r = cellfun(@numel,C.bias);
        if numel(unique(r))>1
            b = nan(max(r),size(C,1));
            s = b;
            for ci = 1:size(b,2)
               s(1:r(ci),ci) = C.speed{ci};
               b(1:r(ci),ci) = C.bias{ci};
            end
            b(b==0) = NaN;
        else
            b = cell2mat(C.bias');
            s = cell2mat(C.speed');
        end
        % exclude bias=0
        s(b==0) = NaN;
        b(b==0) = NaN;
        sm = nanmean(s.*b);
        % for positive values
        if sum(sm>0)>0
            s1 = s;
            s1(b==-1)=NaN;
            C.speedout(sm>0) = nanmean(s1);
            C.biasout(sm>0)=1;
        end
        if any(isinf(C.biasout))
            % for negative values
            if sum(sm<0)>0
                s1 = s;
                s1(b==1)=NaN;
                C.speedout(sm<0) = nanmean(s1);
                C.biasout(sm<0)=-1;
            end
            if any(isinf(C.biasout))
                % for zero values
                if sum(sm==0)>0
                    s1 = s;
                    s1(b~=0)=NaN;
                    C.speedout(sm==0) = nanmean(s1);
                    C.biasout(sm--0)=0;
                end        
            end
        end
        TR.speedout(ind)=C.speedout;
        TR.biasout(ind)=C.biasout;
    end
    
end


%% check if all values entered
if sum(isinf(TR.speedout))~=0 || sum(isinf(TR.biasout))~=0
    error('some values not entered');
end

%% update data
rowval = TR.rowid_out;
rowdelete = cell2mat(TR.rowid);
rowdelete(ismember(rowdelete,rowval)) = [];

%% update data
D.bias(rowval) = TR.biasout; % record bias
D.speed(rowval) = TR.speedout; % record speed
D(rowdelete,:) = []; % delete the rest

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