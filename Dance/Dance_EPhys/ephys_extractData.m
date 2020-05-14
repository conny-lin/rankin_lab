function [Data,MWTDB] = ephys_extractData(pMWT,assaytimes,ISI,preplate,assayTapNumber,varargin)


%% setting


%% GET INFO
legend_gangnam = load('/Users/connylin/Dropbox/Code/Matlab/Library RL/Modules/Chor_output_handling/legend_gangnam.mat');
legend_gangnam = legend_gangnam.legend_gangnam;
msrInterest = {'time','id','bias','speed','tap','midline'};

%% check pMWT files
[~,pMWT_success,~] = getpath2chorfile(pMWT,'Gangnam.mat','reporting',0);
MWTDB = parseMWTinfo(pMWT_success);
pMWT = pMWT_success;
gu = unique(MWTDB.groupname);

%% LOAD AND PROCESS FILES
DataVal = cell(size(pMWT));
% produce habituated response summary per plate
for mwti = 1:numel(pMWT)
    processIntervalReporter(numel(pMWT),1,'MWT',mwti); % reporting
    pmwt = pMWT{mwti}; % get basic info
    % get file path
    pf = char(getpath2chorfile(pMWT(mwti),'Gangnam.mat','reporting',0));
    if exist(pf,'file')~=2; display(pMWT(mwti)); error('no file exist'); end
    % load individual worm data
    Import = load(pf);
    % extract data from relevant times
    DataHab = cell(size(assaytimes,2),1);
    for ti = 1:size(assaytimes,2) % for each tap time
        Data1 = Import.Data(extract_valid_wormtime(Import.time,assaytimes(:,ti)));
        for wrmi = 1:numel(Data1)
            clear D;
            D = convert_import2table(Data1{wrmi},legend_gangnam,'msr',msrInterest);
            D = extract_data_timeofinterest(D.time,D,assaytimes(:,ti));
            % exclude data if all nan
            if sum(isnan(D.bias))==size(D,1); Data1{wrmi} = {}; continue; end; 
            
            % synchronize to tap
            [t,tapN] = syncdata2tap(D.time,D.tap,0);
            if tapN==0
                Data1{wrmi} = {};
            else
                D.time = t;
                Data1{wrmi}  = transform_roundSpeedbytime2(D,'overwriteTime',0);
                if isempty(D)==1; error('flag'); end
            end
        end
        Data1(cellfun(@isempty,Data1)) = []; % clear invalid data
        Data1 = cellfun(@table2array,Data1,'UniformOutput',0); % convert to array
        [r,~] = cellfun(@size,Data1); % find size
        DataHab{ti} = [repmat(ti,sum(r),1) cell2mat(Data1)]; % add time id
    end
    
    %% check consistency of data
    [r,c] = cellfun(@size,DataHab);
    % remove data with no data
    if any(r==0) 
        DataHab(r==0) = [];
    end
    [~,gid] = ismember(MWTDB.groupname(mwti),gu);
    DataVal{mwti} = [repmat(mwti,sum(r),1) repmat(gid,sum(r),1) cell2mat(DataHab)];
end

%% make sure all pMWT contains data
i = cellfun(@isempty,DataVal);
if sum(i) > 0
   warning('the following do not contain data:');
   disp(char(pMWT(i)));
   DataVal(i) = [];
end
DataVal = cell2mat(DataVal);
clear Import DataHab Data Db D;

% convert 2 table
Data = array2table(DataVal,'VariableNames',[{'mwtid','groupid' 'timeid'} msrInterest {'timeround'}]);
clear DataVal D;

% validate data
% create standard times
a= assaytimes(1):.1:assaytimes(2);
a = round((a'-(ISI*(assayTapNumber(1)-1)+preplate)).*10)./10;
a(end) = [];
b = tabulate(Data.timeround);
badtime = b(~ismember(b(:,1),a'),1);
% exclude time
Data(ismember(Data.timeround,badtime),:) = [];

















