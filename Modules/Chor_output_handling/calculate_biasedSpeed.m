function [speed,bias] = calculate_biasedSpeed(D,varargin)


% assaytype = 'poststim'; 
%     need to include first row as pre-tapped response, second row as tapped 
%     response, third+rows as responses for analysis
% assaytype pretim, 
%     include responses before the tap (exclude tap)

%% process varargin
assaytype = 'prestim'; % need to include first row as pre-tapped response, second row as tapped response, third+rows as responses for analysis
% assaytype pretim, include responses before the tap (exclude tap)
vararginProcessor;


%% calculation
speed = Inf; bias = Inf;
switch assaytype
    case 'prestim'
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
        % decide bias and valid speed data
        sval = false(size(D,1),1);
        switch numel(dir)
            case 1 % if all going one direction
                sval(:) = true; 
                b1 = D.bias(1); b1(isnan(b1))=4;
                switch b1
                    case 1; bias=1;
                    case 0; bias = 1;
                    case -1; bias = -1;
                    case 4; bias = NaN;
                end
            case 2 % if two directions
                if ismember(dir,[0 1])
                   bias = 1; sval(:)= true;
                elseif ismember(dir,[0 -1])
                    lastdir = D.bias(end);
                    switch lastdir
                        case 0
                            sval(D.bias==0) = true; bias=1;
                        case -1
                            sval(D.bias==-1) = true; bias=-1;
                        otherwise
                            error('code for this scenario');
                    end
                elseif ismember(dir,[1 -1])
                    lastdir = D.bias(end);
                    switch lastdir
                        case 1
                            sval(D.bias==1) = true; bias=1;
                        case -1
                            sval(D.bias==-1) = true; bias=-1;
                        otherwise
                            error('code for this scenario');
                    end
                else
                    error('code for 2 directions scenario');
                end
            case 3 % if three direction
                sval(isnan(D.bias))=false;
                lastdir = D.bias(end);
                if ismember(dir(2:end),[0 1]) % check last two directions
                    sval(D.bias>=0) = true; bias=1;
                elseif ismember(dir(2:end),[-1 0])
                    switch lastdir
                        case -1; sval(D.bias==-1) = true; bias=-1;
                        case 0; sval(D.bias<=0) = true; bias=-1;
                        otherwise; error('code for this scenario');
                    end
                elseif ismember(dir(2:end),[-1 1])
                    switch lastdir
                        case 1; sval(D.bias==1) = true; bias=1;
                        case -1; sval(D.bias==-1) = true; bias=-1;
                        otherwise; error('code for this scenario');
                    end
                else
                end
            otherwise
                error('more than 3 directions');
        end
        if isnan(bias) && numel(dir)>1; error('deal with mixed NaN dir'); end
        if sum(sval)==0
            speed = NaN;
        else
            speed = nanmean(D.speed(sval)).*bias;    
        end
end


