function MWTSet = Dance_Raster_v1707(pMWT,pSave,varargin)

%% DEFAULT VARIABLES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
graphopt = true; % plot on or off
analysis_purpose = 'rType'; % other options: 'Graph' = generate graphs
displayopt = true;
% tap = 1:30;
% aftertap = the seconds after taps the data analysis program will pull. The smaller this is, the more data will be valid for analysis.
% beforetap = the seconds before taps the data analysis program will pull. The smaller this is, the more data will be valid for analysis.
switch analysis_purpose
    case 'rType' % analyze response types. baseline is determined by movement during 0.6s before a tap and response types are determined by movement within 1s after a tap
        aftertap = 1;
        beforetap = 0.6;
    case 'Graph_Full_Duration' % look at full duration before and after taps
        aftertap = 9;
        beforetap = 1;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% VARARGIN PROCESSOR %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
funname = mfilename;
vararginProcessor
% vararginProcessor_v1707(whos, varargin);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% STANDARD DANCE PROCESSING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[MWTSet,pSave] = DanceM_MWTSetStd_v1707(pMWT,mfilename('fullpath'),varargin,pSave);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% CALCULATE INTPUT DATA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
DB = MWTSet.MWTDB; % get MWT info
% isi
isi = unique(DB.ISI); % get isi information
if numel(isi) ~= 1; error('%s can only accomodate 1 ISI',funname); end
% preplate
preplate = unique(DB.preplate); 
if numel(preplate) ~= 1; error('%s can only accomodate 1 preplate',funname); end
% tapN
tapN = unique(DB.tapN); 
if numel(tapN) ~= 1; error('%s can only accomodate 1 tapN',funname); end
% if tapNmax>tapN; error('# tap requested is more than actual number of taps');end
% taps
taps = 1:tapN;
% create time frames
taptimes = (taps*isi)+preplate-isi;
assayStartTimes = taptimes-beforetap; 
assayEndTimes = taptimes+aftertap;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% CHOR %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[MWTSet,MWTDB] = DanceM_chor(MWTSet,'TrinityOnly'); % chor
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% GRAB DATA  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
gnlist = unique(MWTDB.groupname); % get unique groups 
MasterData = cell(numel(pMWT),numel(gnlist));
fprintf('-- Loading Trinity data, this will take a while ...\n ')
for gi = 1:numel(gnlist) % cycle through groups
    % prepare local variables ---------------------------------------------
    gn = gnlist{gi}; % get group name
    pMWTG = MWTDB.mwtpath(ismember(MWTDB.groupname, gn)); % get pMWT path
    %----------------------------------------------------------------------

    % cycle through plates ------------------------------------------------
    for mwti = 1:numel(pMWTG)
        pMWTp = pMWTG{mwti};
        processIntervalReporter(numel(pMWTG),5,sprintf('group%d plate',gi),mwti)
        TrinityData = load_Trinity_v1707(pMWTp); % load plate data
        D = trinityExtractData_v1707(TrinityData,taps,beforetap,aftertap,{'wormid','time','tap','speed','bias'});
        MasterData{mwti,gi} = D;
    end
    %----------------------------------------------------------------------
end
fprintf('completed\n');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% RASTER PLOT - RED + GREEN + ephys data output
for gi = 1:numel(gnlist) % cycle through groups
    % prepare local variables ---------------------------------------------
    gn = gnlist{gi}; % get group name
    i = ismember(MWTDB.groupname, gn);
    mwtid = MWTDB.mwtid(i);
    MWTDBG = MWTDB(i,:);
    pSaveG = create_savefolder(pSave,gn); % create group folder
    %----------------------------------------------------------------------
    for ti = 1:numel(taps) 

        A = MasterData(:,gi);
        i = cellfun(@isempty,A);
        mwtid = mwtid(~i);
        A(i,:) = [];
        G = cell(size(A,1),1);
        MWTID = G;
        WRMID = G;
        t = nan(size(A,1),2);
        Times = G;
        for mwti = 1:size(A,1)
            D = A{mwti,gi};
            D = D(D.assayperiod == taps(ti),:);

            [T,G{mwti},~] = conv_data_trinity2raster(D,'velocity',{'wormid'});
            [~,Times{mwti},~] = conv_data_trinity2raster(D,'time',{'wormid'});

            WRMID{mwti} = T.wormid;
            MWTID{mwti} = repmat(mwtid(mwti),size(Array,1),1);
            
            t(mwti,1) = min(D.time);
            t(mwti,2) = max(D.time);
        end
        Data = cell2mat(G);
        mwtid = cell2mat(WRMID);
        rTime = cell2mat(Times);
        
        % save data
        n = size(D,1);
        t1 = min(t(:,1));
        tf = max(t(:,2));
        savename = sprintf('raster_%.0fs_%.0fs_N%d',t1,tf,n);
        p = fullfile(pSaveG, sprintf('%s.mat',savename));
        save(p,'MWTDBG','Data','mwtid','rTime');
        % ----------------------------------


        %% graphs ---------------------------------------------------------
        if graphopt
            % make hot figure
            f1 = rasterPlot_colorSpeed_hot(Data);
            p = create_savefolder(fullfile(pSaveG,'raster hot'));  % get print path
            set(f1,'PaperPositionMode','auto'); % set to save as appeared on screen
            cd(p); print (f1,'-depsc', '-r1200', savename); % save as eps
            print(f1,'-dtiff', savename); % save as (r600 will give better resolution)
            close;


            % make figure (green, jet)
            [f1,~] = rasterPlot_colorSpeed_gradient2(Data,rTime,0);
            p = create_savefolder(fullfile(pSaveG,'raster jet'));  % get print path
            % print
            set(f1,'PaperPositionMode','auto'); % set to save as appeared on screen
            cd(p); print (f1,'-depsc', '-r1200', savename); % save as eps
            print(f1,'-dtiff', savename); % save as (r600 will give better resolution)
            close;


            % make lava plot (sort proportion)
            f1 = moveproportion(Data,rTime,0);
            p = create_savefolder(fullfile(pSaveG,'lava proportion plot'));  % get print path
            titlename = sprintf('%s/speed_dir_proportion_%s',p,savename);
            print(f1,'-depsc', titlename); 
            print(f1,'-dtiff', titlename); % save as (r600 will give better resolution)
            close;
        end
    end

end



