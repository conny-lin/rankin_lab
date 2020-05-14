function [A,N,NV,plateidlist] = responseType_pct_plate(plateinfo,R,responseType,varargin)

%% revision
% r20170519:
% - removed worms that has no data across all time assayed
% - fill in no data as empty cell


%%
Type = 'per time'; % or per time
vararginProcessor
%%


%% SUMMARIZE RESPONSE TYPE PER TIME PER GROUP PER PLATE
plateidu = unique(plateinfo.mwtid);
switch Type
    case 'per time'
        A = nan(size(plateidu,1), size(R,2));
    case 'any'
        A = nan(size(plateidu,1), 1);
end

N = A;
NV = A;
plateidlist = nan(size(plateidu,1),1);


for plateidi = 1:numel(plateidu)
    % get plate id
    pid = plateidu(plateidi);
    i_plate = plateinfo.mwtid == pid;
    
    % get data for this plate
    R1 = R(i_plate,:);
    R1 = table2cell(R1);
    R1(sum(cellfun(@isempty,R1),2) == size(R1,2),:) = []; % r20170519: exclude worms that does not have data for all duration
    %     R1 = regexprep(R1,'','no data'); % r20170519
    R1(cellfun(@isempty,R1)) = {''}; % r20170519: fill in no data as empty vector
%%
    
    switch Type
        case 'per time'
             R1 = categorical(R1);
             if size(R1,1) == 1 % if there is only one plate for this group
                p = ones(1,5);
                n = p;
                n2 = n;

            else

                % get counts per category
                d = countcats(R1);
                leg = categories(R1);

                r = d(ismember(leg, responseType),:); % define response type
                if isempty(r)
                   r = zeros(1, size(R1,2));
                end

                n = sum(d,1);

                r_exclude = d(ismember(leg, 'baseline'),:);
                if ~isempty(r_exclude)
                    n2 = (n - r_exclude);
                else
                    n2 = n;
                end
                p =  r./n2;
             end
            
             
        case 'any'
            
            n2 = sum(any(ismember(R1,responseType),2)); % count any time the worm response as responseType

            n = size(R1,1);
            p = n2/n;
    end
    
    
    % summarize
    A(plateidi,:) = p;
    NV(plateidi,:) = n2;
    N(plateidi,:) = n;
    plateidlist(plateidi) = pid;
end