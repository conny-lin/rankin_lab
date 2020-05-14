function MWTSet = Dance_Rave(pMWT,varargin)
%% Dance_Glee_Acafellas     make raster plot for pMWTS by group
% INPUT OPTIONS:
%     pSave
%     analysisNLimit = 5; put minimum number of worms must be reversing to
%     first tap in order to qualify for analysis
%     alcoholTestDose
%     refStrain = 'N2';



%% DEFAULTS 
nInput = 1;
pData = '/Volumes/COBOLT/MWT';
pSave = '/Users/connylin/Dropbox/RL/Dance Output';
chorName = 'Trinity';
pAnalysis = pData;
% set analysis time seeting
startTime = 88;
endTime = 405;
assayInt = 10;
assaywindow = 10;
frameInt = 0.2;
% analysis restrictions
analysisNLimit = 5; % number of worms must exist before including in analysis
alcoholTestDose = 400;
refStrain = 'N2';
tapEqual = true; % include only groups with equal tap number use 0 to ignore
NWORMS = Inf;
% stats setttings
alphaValue = 0.05;
posttestname = 'bonferroni';
% graphic settings
visibleG = 0;
% data processing settings
saveoption = 0; % do not save trinity output
cleanup = 1;
alignTap = 1;
%% STANDARD DANCE PROCESSING
% add general function paths
pDanceM = [fileparts(fileparts(mfilename('fullpath'))),'/Modules'];
addpath(pDanceM); DanceM_funpath; 
% VARARGIN PROCESSOR
if nargin > nInput
    A = varargin;
    if isinteger(numel(A)/2) == 1
        error('variable entries incorrect')
    end
    callName = A(1:2:numel(A));
    % check if call names are strings
    for x = 1:numel(callName)
       if ischar(callName{x}) == 0
           error('varargin incorrect');
       end
    end
    % get values
    setValue = A(2:2:numel(A));
    for x = 1:numel(callName)
        % determine input type
        if eval(sprintf('ischar(%s)',callName{x})) == 1
            eval(sprintf('%s = ''%s'';',callName{x},setValue{x}));
        elseif eval(sprintf('isnumeric(%s)',callName{x}))
            eval(sprintf('%s = %d;',callName{x},cell2mat(setValue(x))));
        elseif eval(sprintf('iscell(%s)',callName{x}))
            a = setValue{x};
            eval(sprintf('%s = a;',callName{x}));
        else
            error('need coding');
        end
    end
end
% create MWTSet
MWTSet.ProgramName = mfilename;
MWTSet.timestamp = generatetimestamp; % generate timestamp
MWTSet.PATHS.pMWT = pMWT;
% save folder
pSave = [pSave,'/',mfilename];
if isdir(pSave) == 0; mkdir(pSave); end
MWTSet.PATHS.pSaveA = pSave;
MWTSet = DanceM_createIndex(MWTSet); % create index


return
%% DANCE STANDARD PROCSSING -----------------------------------------------
%% EVALUATE MWT INPUT
% if tapEqual == true && numel(MWTSet.Info.VarIndex.tapN) > 1
%     error('tapN is not equal across inputs, set tapEqual = 0 to over come this');
% elseif tapEqual == true && numel(MWTSet.Info.VarIndex.tapN) == 1
%     tapNExpected = unique(MWTSet.Info.VarIndex.tapN);
% end
%% CHOR
pMWTc = convertTrinityDat2Mat(pMWT,1); 
L = chormaster4(chorName,pMWTc);
% summarize trinity data and delete .dat file to save memory
pMWTc = convertTrinityDat2Mat(pMWTc,1); 


%% DATA LOADING -----------------------------------------------------------
%% IMPORT RAW DATA & COMPILE TO .MAT (.trinity.*.dat)
% find and load trinity summary file, optional clean up
MDb = MWTSet.Info.MWTDb;
pMWT = MDb.mwtpath;
fMWT = MDb.mwtname;
Imp = table;
Imp.mwtid = MWTSet.Info.MWTDbInd.mwtid;
Imp.data = cell(size(MDb,1),1);
for mwti = 1:numel(pMWT)
    pMWTp = pMWT{mwti};
    if isempty(dircontent(pMWTp,'trinitySummary.mat')) == 1
       % if can't find trinitySummary.mat file, deal with errors
       fprintf('no trinitySummary.mat file found in [%s]\n',fMWT);       
    else
        % load file
        D = load([pMWTp,'/trinitySummary.mat'],'masterData'); D = D.masterData;
        Imp.data{mwti} = D;
    end
end
% record problem import
i = cellfun(@isempty,Imp.data);
MWTSet.Info.MWTDb.Import_trinity = i;
% delete problem import
Imp(i,:) = [];


% legend
importname = 'trinity';
Legend = {'time' 'number' 'goodnumber' 'speed' 'speed_std'...
            'bias' 'tap' 'puff' 'loc_x' 'loc_y' ...
            'morphwidth' 'midline' 'area' 'angular' 'aspect'...
            'kink' 'curve' 'crab'};
% store in MWTSet
MWTSet.Info.MsrList.(importname) = Legend';


% translate get column 1(time), 4(speed), 6(bias)
% translate import into tables
% create table variablenames
a = Imp.Properties.VariableNames;
L = [a(1:end-1) {'wormid'}  Legend];
nCol = numel(L);
T = table;
for mwti = 1:size(Imp,1)
    % get data
    D = Imp.data{mwti}; 
    nworm = size(D,1); % get number of worms
    wormidList = cellfun(@str2num,D(:,1)); % get wormid
    % get MWTindex
    ind = table2array(Imp(mwti,1:end-1));
    % create output array
    [r,~] = cellfun(@size,D(:,2));
    nRow = sum(r);
    A = nan(nRow,nCol);
    r1 = 1;
    for wi = 1:nworm
        r2 = r1+r(wi)-1; % get end row 
        d = D{wi,2}; % get data
        wormid = repmat(wormidList(wi),size(d,1),1); % get worm id
        v = repmat(ind,size(d,1),1); % variable index
        d = [v wormid d]; % add to d
        A(r1:r2,1:nCol) = d; % store in A
        r1 = r2+1; % reset start row
    end
    T = [T;array2table(A,'VariableNames',L)]; % store in master table;
end
%% remove invalid data
T(isnan(T.speed) | isnan(T.bias) |T.time < startTime | T.time > endTime,:) = [];
MWTData.trinity = T;
clearvars Imp T A D d; % clean memory




%% CALCULATE ANALYSIS SETTING BASED ON DATA LOADED ------------------------
%% ASSAY TIMES
D = MWTData.trinity;
%% create standard assay time based on setting
startList = startTime:assayInt:endTime;
endList = startList+assaywindow;
assayBlocksStd = [startList;endList]; 


% create assay times  - if no tap align
if alignTap == 0 || sum(D.tap) == 0 %% align to tap if required to do so
    T = nan(size(D,1),1);
    time = D.time;
    for ti = 1:size(assayBlocksStd,2)
        i = time >= assayBlocksStd(1,ti) & time < assayBlocksStd(2,ti);
        T(i) = assayBlocksStd(1,ti);
    end
    D.assayTime = T;
end

% create assay times  - if tap align
if alignTap == 1 && sum(D.tap) > 0 % if there is tap
    fprintf('Translating time to assay times, this will take a while\n');
    % detect if there is tap
    mwtu = unique(D.mwtid);
    D.assaytimeseq = nan(size(D,1),1);
    D.assayframetime = nan(size(D,1),1);
    D.assayblock = nan(size(D,1),1);
    for mwti = 1:numel(mwtu) % for each plate
        fprintf('-%d/%d mwt',mwti,numel(mwtu));
        % find tap time
        tapTime = unique(D.time(D.tap == 1 & D.mwtid == mwtu(mwti)));
        imwt = D.mwtid == mwtu(mwti); % index to mwt plate
        for tbi = 2%1:size(assayBlocksStd,2) % for each assay block
            start = assayBlocksStd(1,tbi);
            finish = assayBlocksStd(2,tbi);
            i = tapTime >= start & tapTime < finish;
            if isempty(t) == 0 % if one or more tap is found, 
                tstart = t(1);  % align to the earlier tap
                assayTime = [flip(tstart-frameInt:-frameInt:start)  tstart:frameInt:finish+frameInt];
            elseif sum(i) == 0 % if no tap found in the time frame, 
                tstart = start + tapTime(1) - floor(tapTime(1)); % align to first tap
                assayTime = [flip(tstart-frameInt:-frameInt:start)  tstart:frameInt:finish+frameInt];
            end
            %% translate frametime(ft) into assay times
            j = imwt & D.time >= start & D.time < finish;
            d = D.time(j);
            s = repmat(assayTime(1:end-1),size(d,1),1);
            f = repmat(assayTime(2:end),size(d,1),1);
            t = repmat(d,1,numel(assayTime)-1);
            v = t >= s & t < f; % later than 
            ts = nan(size(d,1),1);
            ft = nan(size(d,1),1);
            ab = nan(size(d,1),1);
            for x = 1:size(v,1)
                a = find(v(x,:));
                if isempty(a) == 0
                   ts(x) = a;
                   ft(x) = assayTime(a);
                   ab(x) = tbi;
                elseif numel(a) > 1
                    error('too many qualified');
                end
            end
            D.assaytimeseq(j) = ts;
            D.assayframetime(j) = tf;
            D.assayblock(j) = tb;
        end
    end
end
%                 
%             
%             return
%         
%             if isempty(i) == 0 % if tap is found
%                 % find taps within start and finish time
%                 t = tapTime(:,mwti);
%                 i = t > start & t < finish;
%                 t = t(i);
%                 if isempty(t) == 0 % if one or more tap is found, 
%                     % align to the first tap
%                     tstart = t(1);
%                     assayTime = [flip(tstart-frameInt:-frameInt:start)  tstart:frameInt:finish+frameInt];
%                 elseif sum(i) == 0 % if no tap found in the time frame, 
%                     % align to first tap
%                     df = i(1) - floor(i(1));
%                     tstart = start + df;
%                     assayTime = [flip(tstart-frameInt:-frameInt:start)  tstart:frameInt:finish+frameInt];
%                 end
%             else % if no tap is found
%                 assayTime = start:frameInt:finish;
%             end
%             assayFrames{mwti,ti} = assayTime;   
%         end
%         return
% %         % separate tap time for each plate
% %         mwtu = unique(u(:,1));
% %         A = [];
% %         a = u(u(:,1) == mwtu(x),2);
% %         if x == 1
% %             A = a;
% %         else
% %             if numel(a) == size(A,1)
% %                 A = [A a];
% %             else
% %                 error('tap time inconsistent');
% %             end
% %         end
%     end
%     
% 
%     %% make assay time based on aligned tap
%     assayFrames = cell(numel(unique(D.mwtid)),size(assayBlocksStd,2));
%     mwtu = unique(D.mwtid);
%     for mwti = 1:numel(mwtu)
%         for ti = 1:size(assayBlocksStd,2)
%             start = assayBlocksStd(1,ti);
%             finish = assayBlocksStd(2,ti);
%             if isempty(tapTime) == 0 % if tap is found
%                 % find taps within start and finish time
%                 t = tapTime(:,mwti);
%                 i = t > start & t < finish;
%                 t = t(i);
%                 if isempty(t) == 0 % if one or more tap is found, 
%                     % align to the first tap
%                     tstart = t(1);
%                     assayTime = [flip(tstart-frameInt:-frameInt:start)  tstart:frameInt:finish+frameInt];
%                 elseif sum(i) == 0 % if no tap found in the time frame, 
%                     % align to first tap
%                     df = tapTime(1) - floor(tapTime(1));
%                     tstart = start + df;
%                     assayTime = [flip(tstart-frameInt:-frameInt:start)  tstart:frameInt:finish+frameInt];
%                 end
%             else % if no tap is found
%                 assayTime = start:frameInt:finish;
%             end
%             assayFrames{mwti,ti} = assayTime;   
%         end
%     end
%     assayFrames_mwtname =  mwtu;




%% VALIDATE DATA
fprintf('validating data, this will take a while\n');
% assayTime must has a value
D(isnan(D.assaytime),:) = [];
% bias can not be nan
% exclude worms not exist for fll duration within time frame
a = unique([D.mwtid D.wormid D.assayblock],'rows');
ID = array2table(a,'VariableNames',{'mwtid','wormid','assayblock'});
val = true;
nTimeVal = numel(unique(D.assaytimeseq));
reportInt = 1:100:size(ID,1);
for id = 1:size(ID,1)
    if sum(ismember(id,reportInt)) == 1
        fprintf('-check set %d/%d\n',id,size(ID,1));
    end
    block = ID.assayblock(id);
    mwtid = ID.mwtid(id);
    r = assayFrames_mwtname == mwtid;
    i = D.mwtid == ID.mwtid(id) & D.wormid == ID.wormid(id) &...
        D.assayblock == block;
    if numel(unique(D.assaytimeseq(i))) ~= nTimeVal
        val(i) = false; % mark as bad
    end
end
return

%%
D(~val,:) = [];
MWTData.trinity = D;
fprintf('\n');


%% CALCULATION: RASTER PLOT BIAS* SPEED

%% average speed * bias per frametimes(ftU)
% PROBLEM: sometimes speed has value as large as 1.32 but bias is 0. What
% does that mean?

% take data needed
D = MWTData.trinity
% exclude bias == nan
D(isnan(D.bias),:) = [];
% calculate speed dir
D.speedDir = D.speed.*D.bias;  
%%
a = unique([D.groupname D.wormid D.assayblock],'rows');
ID = array2table(a,'VariableNames',{'mwtid','wormid','assayblock'});
gnu = unique(D.groupname);
for gi = 1:numel(gnu)
%    for  
end

%%

for x = 1:numel(ftU)
    t = ftU(x);

    d = D(D.frametime == t & D.wormid == wormidU(wi),:);
    if sum(diff(d.bias)) == 0 
    % if move within this frame all in one dir
    % calculate mean of all data
        A(wi,x) = mean(d.speedDir);
    elseif sum(diff(d.bias)) ~= 0 
    % if movement within this frame contains dir shift
        if sum(d.bias < 0) > 0
        % if reversal exists, only calculate reversal mean
            A(wi,x) = mean(d.speedDir(d.bias < 0));
        elseif sum(d.bias > 0) > 0
        % if only forward exists, only calculate forward mean
            A(wi,x) = mean(d.speedDir(d.bias > 0));
        elseif sum(d.bias == 0) > 0
        % if only bias == 0 (presumably pause), only calculate pase
        % speedDir, which is set to zero
            A(wi,x) = mean(d.speedDir(d.bias == 0));
        end
    end
end
plateSumm = A;
% if any nan results, flag
if sum(any(isnan(A))) > 0
   error('some nan results, code to fix'); 
end

% text output
if saveoption == 1
    dlmwrite(sprintf('%s/%ds_%ds_fint_%.1f_N%d_wormID.trinitySummary',pMWTp,start,finish,frameInt,size(wormID,1)),wormID);
    dlmwrite(sprintf('%s/%ds_%ds_fint_%.1f_N%d_frametime.trinitySummary',pMWTp,start,finish,frameInt,size(wormID,1)),frameTimeRecord)
    dlmwrite(sprintf('%s/%ds_%ds_fint_%.1f_N%d_rasterdata.trinitySummary',pMWTp,start,finish,frameInt,size(wormID,1)),plateSumm)
end
    

% %% survey times
% %  PROBLEM: bias can be -1, 1, or 0. 
% % check if nan for bias
% D = DataAll;
% if sum(isnan(D.bias)) > 0
%     error('some bias record are NaN');
% end
% % delete data outside of assay times
% D(D.time < assayTime(1) | D.time >= assayTime(end),:) = []; 
% 
% 
% % translate frametime(ft) into assay times
% ft = D.time;
% for ti = 1:numel(assayTime)-1
% %     if ti == numel(assayTime)
% %         i = ft >= assayTime(ti);
% %     else
%         i = ft >= assayTime(ti) & ft < assayTime(ti+1);
% %     end
%     if sum(i) == 0
%        [~,fn] = fileparts(pMWTp);
%        plateSumm = [];
%        warning('time point (%.2fs) in plate [%s] has no data, skip',assayTime(ti),fn); 
%        return
%     end
%     ft(i) = assayTime(ti);
% end
% % validate:
% % 1. all frame time must equal to one of the assay time
% % 2. number of unique frame time must equal to (assay time -1)
% if sum(ismember(ft,assayTime)) ~= numel(ft) ||...
%         numel(unique(ft)) ~= numel(assayTime)-1
%     error('frame time matching failed');
% else
%     D.frametime = ft;
% end
% 
% % all worms must exist throughout the period
% widU = unique(D.wormid);
% t1 = assayTime(1);
% t2 = assayTime(end-1);
% n = nan(numel(widU),1);
% nV = (numel(assayTime)-1);
% for wi = 1:numel(widU)
%     i = D.wormid == widU(wi);
%     t = D.frametime(i);
%     n = numel(unique(t));
%     if n~=nV
%         D(i,:) = []; % delete
%     end
% end



%% GRAPH: RASTER
% run raster plot
for gi = 1:numel(groupnameList)
    gT = groupnameList{gi};
    pMWT = DbT.mwtpath(ismember(DbT.groupname,gT));
    % raster plot per time points
    for ti = 1:numel(startList)
        % create save folder
        pSave = sprintf('%s/%s/%s/%s',pDest,strainName,gT);
        if isdir(pSave) == 0; mkdir(pSave); end
        % run raster
        start = startList(ti);
        finish = endList(ti);
        [f1,Data,savename,Tog,mwtid,rTime] = rasterPlot_colorSpeed(pMWT,start,finish,...
            'NWORMS',Inf,'visibleG',0,'frameInt',frameInt);
        % save fig
        cd(pSave);
        set(f1,'PaperPositionMode','auto'); % set to save as appeared on screen
        print (f1,'-depsc', '-r1200', savename); % save as eps
        close;
        % save data
        save(sprintf('%s/%s rasterData.mat',pSave,savename),'Data','mwtid','Tog','rTime');
    end
end

% summarizing
% forward movement
For = plateSumm;
For(For<0)=0; % zeros for reversal
For = For + 0.2; % add 0.2
For(For == 0.2)= 0;

% backwards
Bac = plateSumm;
Bac(Bac>0)=0;
Bac = Bac - 0.4;
Bac(Bac==(-0.4))=0;

% Tog
Tog = For + Bac;
Tog(Tog>0.8)=0.8; % max out at 0.8
Tog(Tog<-0.8)=-0.8; % max out at -0.8
Tog(Tog==0)=0.2; % if zero then 0.2
Tog(1,1)=0.8; 
Tog(1,2)=-0.8;
 

% scale output to N = scaleNumWorms
nWorms = size(Tog,1);
if isinf(NWORMS) == 0
    if nWorms < NWORMS % if less worm than defined in scaleNumWorms, add zeros
        extr = NWORMS - nWorms;
        extrRows = zeros(extr, size(Tog,2));
        Tog = [extrRows;Tog]; 
    else % take the first 200 worms
        Tog(NWORMS:end,:)=[];
    end
end



% create image
if visibleG == 0
    figure1 = figure('Color',[1 1 1],'Visible','off');
else
    figure1 = figure('Color',[1 1 1],'Visible','on');
end
axes1 = axes('Parent',figure1,...
    'ZColor',[0.8 0.8 0.8],...
    'YDir','reverse',...
    'YColor',[0.8 0.8 0.8],...
    'XColor',[0.8 0.8 0.8],...
    'Layer','top');
imagesc(Tog,'Parent',axes1)


% save
% create save name
starttxt = num2str(start);
finishtxt =  num2str(round(finish));
ntext = num2str(size(Tog,1));
savename = ['rasterPlot_',starttxt,'_',finishtxt,'_N_',ntext];

%% ORGANIZE OUTPUT BY GROUP
% GET GROUP NAME AND GRAPHING SEQUENCE
% p = MWTSet.Data.Raw.pMWT;
% [~,GroupName] = cellfun(@fileparts,cellfun(@fileparts,p,'UniformOutput',0),'UniformOutput',0);
% gnameU = unique(GroupName);
% i = ~cellfun(@isempty,regexp(gnameU(:,1),'^N2')); % N2 graph first
% GroupSeq = gnameU([find(i);find(~i)]);
% 
% [~,fMWT] = cellfun(@fileparts,p,'UniformOutput',0);
% D = MWTSet.Data.Raw;
% A = struct;
% for g = 1:size(GroupSeq,1)
%     B = struct;
%     gname = GroupSeq{g};
%     i = ismember(GroupName,GroupSeq(g));
%     B.MWTplateID = fMWT(i);
%     B.MWTind = find(i);
%     B.time = D.X(:,i);
%     B.N_NoResponse = D.N.NoResponse(:,i);
%     B.N_Reversed = D.N.Reversed(:,i);
%     B.N_TotalN = D.N.TotalN(:,i);
%     B.RevFreq_Mean = D.Y.RevFreq(:,i);
%     B.RevFreq_SE = D.E.RevFreq(:,i);
%     B.RevDur_Mean = D.Y.RevDur(:,i);
%     B.RevDur_SE = D.E.RevDur(:,i);
%     B.RevSpeed_Mean = D.Y.RevSpeed(:,i);
%     B.RevSpeed_SE = D.E.RevSpeed(:,i);
%     A.(gname) = B;
% end
% MWTSet.Data.ByGroupPerPlate = A;



%% STANDARD INFORMATION OUTPUT
% expinfosheet
writetable(MWTSet.Info.MWTDb,'expinfosheet.csv');

% mat files
cd(pSave); save(mfilename,'MWTSet');


%% Report done
fprintf('\n\n***DONE***\n\n');
















