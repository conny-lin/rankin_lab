function MWTSet = Dance_rType3_v1707(pMWT,pSave,varargin)
%% notes
% baseline assay 0.3s before tap
% response assayed till 1s after tap
% above defined in getAssayTime
% frameint = 0.05;

%% DEFAULT VARIABLES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% standard dance vairables
displayopt = false;
overwrite = false;
timename = 'timed';
varnames = {'speed','bias'};
frameint = 0.05;
msrlist = {'F','R','P'};
rTargetName = {'acc','rev','pause'};
aggregateTypes = {'F',{'forward','forward_acc'};
    'R',{'reverse','reverse_acc'};
    'P',{'pause','forward_dec','reverse_dec'}};
% type = 'tap';
% n_lowest = 10;

% rTarget = {'accelerate forward'; 'reversal';'pause'};
genotype = '';
strain = '';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% VARARGIN PROCESSOR %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
vararginProcessor;  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% STANDARD DANCE PROCESSING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[MWTSet,pSave] = DanceM_MWTSetStd_v1707(pMWT,mfilename('fullpath'),varargin,pSave);
[MWTSet,MWTDB] = DanceM_chor(MWTSet,'TrinityOnly','displayopt',displayopt,'overwrite',overwrite); % chor
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% get information from trinity files %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
TrinityInfo = getinfo_trinityind(MWTDB);
TimeSet = getAssayTime(MWTDB,'rType'); % get time assay information
[TrinityInfo,TimeSet] = get_trinityIndInfo(TrinityInfo,TimeSet);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% get data from trinity %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
gnlist = unique(MWTDB.groupname);
tN =numel(TimeSet.atstart);
OUT = struct;
fnames = {'mwtid','wormid','groupname','time','rtype'};
for fi =1:numel(fnames)
    OUT.(fnames{fi}) = cell(tN,1);
end
% DATA = struct;

outputrow = 1;

for gi = 1:numel(gnlist) % cycle through groups        
    gn = gnlist{gi};
    for tii = 1:tN % cycle through time points 
        
        loopreporter(tii,sprintf('group%d time',gi),1,tN);

        % get time ----------------------------------------------------
        starttime = TimeSet.atstart(tii); 
        finishtime = TimeSet.atend(tii); 
        %--------------------------------------------------------------

        % load trinity data -------------------------------------------
        Tinfo = TrinityInfo(TrinityInfo.t1 <= starttime & ...
            TrinityInfo.t2 >= finishtime & ...
            ismember(TrinityInfo.groupname, gnlist{gi}),...
            {'mwtid','mwtpath','filename'});
        pTri = cellfun(@fullfile, Tinfo.mwtpath, Tinfo.filename,'UniformOutput',0); % get trinity file paths
        Tri = load_trinityInd_v1707(pTri, Tinfo.mwtid); % load
        % delete data outside of assay times 
        Tri = trim_trinitydata(Tri,starttime,finishtime,...
            {'id','ids','mwtid','wormid','time','speed','bias','tap'},...
            'removenan',true);
        %--------------------------------------------------------------
        

        %% treat trinity data  -----------------------------------------
        Tri = taptime_align_timed(Tri,timename);
        S = tsfTri2timetable(Tri,varnames,timename,frameint); % create timetable
        S = deal_dirshift(S,timename);
        Data.(gn)(tii) = S;
%         [S,~] = cal_velocity(S);
        % -------------------------------------------------------------

        % get response type  -----------------------------------------
        A = S.speed.timetable; % get data out
        Baseline_speed = table2array(A(A.Time<0,:)).'; % get baseline
        Response_speed = table2array(A(A.Time>0,:)).'; % get response
        A = S.bias.timetable; % get data out
        Baseline_bias = table2array(A(A.Time<0,:)).'; % get baseline
        Response_bias = table2array(A(A.Time>0,:)).'; % get response
        [~,rDecisionText,~] = compute_response_type_v1707(Baseline_speed,Baseline_bias, Response_speed, Response_bias);    
        % -------------------------------------------------------------
        
        % extract plate information
        A = S.speed.timetable;
        a = A.Properties.VariableNames;
        a = regexpcellout(a,'_','split');

        n = size(rDecisionText,1);
        
        % compile data in output structure
        OUT.mwtid{outputrow} = cellfun(@str2num,a(:,2));
        OUT.wormid{outputrow} = cellfun(@str2num,a(:,3));
        OUT.rtype{outputrow} = rDecisionText;
%         OUT.groupname{outputrow} = repmat({gn},n,1);
        OUT.time{outputrow}= repmat(tii,n,1);

        outputrow = outputrow+1;


    end
end

% store
MWTSet.Raw = Data;

% compile 
T = table;
T.mwtid = cell2mat(OUT.mwtid);
T.wormid = cell2mat(OUT.wormid);
% T.groupname = celltakeout(OUT.groupname);
T.time = cell2mat(OUT.time);
T.rtype = celltakeout(OUT.rtype);

MWTSet.Data = T;
save(fullfile(pSave,'data.mat'),'MWTSet');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% PREPARE STATS TABLE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get data
D = MWTSet.Data;
M = MWTSet.MWTDB(:,{'mwtid','groupname'});
D = innerjoin(M,D);
[~,~,Legend] = compute_response_type_v1707(Baseline_speed,Baseline_bias, Response_speed, Response_bias);    

% summarize
timelist = unique(D.time);
gnlist = unique(D.groupname);
T = table;
T1 = table;
for gi = 1:numel(gnlist) % cycle through groups        
    gn = gnlist{gi};
    for tii = 1:numel(timelist) % cycle through time points 
        D1 = D(ismember(D.groupname,gn) & ismember(D.time,timelist(tii)),:);

        % make perecentage table
        [A,~,~,l] = crosstab(D1.mwtid, D1.rtype);
        rnames = l(:,2);
        % create standard output
        A1 = zeros(size(A,1),numel(Legend.name));
        [i,j] = ismember(Legend.name,rnames);
        A1(:,i) = A(:,j(i));  % rearrange
        P = A1./sum(A1,2);
        mwtid = l(:,1);
        mwtid(cellfun(@isempty,mwtid)) = [];
        mwtid = cellfun(@str2num,mwtid);
        t = repmat(tii,numel(mwtid),1);
        A = array2table([mwtid t P],'VariableNames',[{'mwtid','time'} Legend.name]);        
        T = [T;A];
    end
end

MWTSet.DataT = T;
save(fullfile(pSave,'data.mat'),'MWTSet');


%% SUMMARIZE
M = MWTSet.MWTDB(:,{'mwtid','groupname','strain','rx'});
D = innerjoin(M,T);
rx = M.rx;
ictrl = strcmp(rx,'NA');
rxe = rx(~ictrl);
if sum(regexpcellout(rxe,'mM')) == numel(rxe)
    M.dose = rx;
    M.dose(ismember(M.dose,'NA')) = {'0mM'};
    M.rx = [];
else
    M.rx = rx;
end
% examine variable types
factors = {};
if ismember('dose',M.Properties.VariableNames)
    factors{end+1} = 'dose';
else
    factors{end+1} = 'rx';
end
if numel(unique(M.strain)) >1
    factors{end+1} = 'strain';
end

% join table
D = innerjoin(M,D);


%% AGGREGATE TYPES

for msri = 1:size(aggregateTypes,1)
    msr = aggregateTypes{msri,1};
    a = [];
    for ti = 1:numel(aggregateTypes{msri,2})
        a = a+aggregateTypes{msri,2}{ti};
    end
    D.(msr) = a;
end


MWTSet.Data_PCT = D;
save(fullfile(pSave,'data.mat'),'MWTSet');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% DESCRIPTIVE STATS & RMANOVA N=PLATES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% anova settings +++++++++++++++++
MWTDB = MWTSet.MWTDB;
rmName = 'time';
idname = 'mwtid';
Data = MWTSet.Data_PCT;

% descriptive stats +++++++++++++++++++++++++++++++
txt1 = '*** Descriptive Stats ***';
t = tabulate(MWTDB.groupname);
% gname
gname = sortN2first(t(:,1),t(:,1));
gname = strjoin(gname,', ');
txt2 = sprintf('gname = %s',gname);
% n
n = sortN2first(t(:,1),t(:,2));
n = strjoin(cellfun(@num2cellstr,n),', ');
txt3 = sprintf('N(plates) = %s',n);
Ntext = sprintf('%s\n%s\n%s',txt1,txt2,txt3);
% -----------------------------------------------

for msri = 1:numel(msrlist)

    msr = msrlist{msri};

    % anova +++++++++++++++++++++++++++++++++++++++
%     A = anovarm_convtData2Input(Data,rmName,[fName,fNameInd],msr,idname);
%     textanova = anovarm_std(A,rmName,fNameInd,fNameS);
  [textanova,DS] = anovarm_std(Data,rmName,factors,idname,msr);
    % -----------------------------------------------

    textout = sprintf('%s\n%s',Ntext,textanova); % join
    
    % save data ++++++++++++
    cd(pSave);
    p = sprintf('%s RMANOVA.txt',msr);
    fid = fopen(p,'w');
    fprintf(fid,'%s',textout);
    fclose(fid);
    % ----------------------

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% graph %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
G = grpstats(Data,{'groupname','time'},{'numel','mean','sem'},'DataVars',msrlist);
MWTSet.Descriptive = G;
% store in struct
S = struct;
gnlist = unique(G.groupname);
gnlist = sortN2first(gnlist,gnlist);

for msri = 1:numel(msrlist)
    msr = msrlist{msri};
    A = nan(numel(unique(G.time)), numel(gnlist));
    S.(msr).gn = gnlist;
    S.(msr).N = A;
    S.(msr).X = A;
    S.(msr).Y = A;
    S.(msr).E = A;

    for gi = 1:numel(gnlist)
        gn = gnlist{gi};
        D = G(ismember(G.groupname, gn),:);
        r = D.time;
        S.(msr).Y(r,gi) = D.(sprintf('mean_%s',msr));
        S.(msr).E(r,gi) = D.(sprintf('sem_%s',msr));
        S.(msr).N(r,gi) = D.(sprintf('numel_%s',msr));
        S.(msr).X(r,gi) = r;
    end
end
MWTSet.Graph = S;

%% Make graphs
G = MWTSet.Graph;
for msri = 1:numel(msrlist)

    msr = msrlist{msri};
    
    % get graph data ----------------
    D = G.(msr);
    Y = D.Y;
    E = D.E;
    X = D.X;
    gn = D.gn;
    % -------------------------------

    % setting --------------------------
    w = 4.5;
    h = 3.5;
    gp = graphsetpack('cathyline');
    gnss = regexprep(gn','_',' ');
    gnss = regexprep(gnss,strain,genotype);
    gnss = regexprep(gnss,'N2','Wildtype');
    gp.DisplayName = gnss;
    gpn = fieldnames(gp);
    % -----------------------------------

    % plot ------------------------------------------
    fig = figure('Visible','off','PaperSize',[w h],'Unit','Inches',...
        'Position',[0 0 w h]); 
    ax1 = axes('Parent',fig,'XTick',[0:10:30]); 
    hold(ax1,'on');
    e1 = errorbar(X,Y,E,'Marker','o','MarkerSize',4.5);
    % -----------------------------------------------

    % get settings -----------------------------------
    for gi = 1:numel(gn)
        for gpi = 1:numel(gpn)
            e1(gi).(gpn{gpi}) = gp.(gpn{gpi}){gi};
        end
    end
    % -----------------------------------------

    % title -------------------------
    title(genotype)
    % -------------------------------

    % legend ------------------------
    %     lg = legend('show');
    %     lg.Location = 'eastoutside';
    %     lg.Box = 'off';
    %     lg.Color = 'none';
    %     lg.EdgeColor = 'none';
    %     lg.FontSize = 8;
    % -------------------------------

    % x axis -----------------------------------------------
    xlabel('Tap')
    xlim([0 30.5]);
    % --------------------------
    % y axis ---------------------
    yname = sprintf('P (%s)',rTargetName{msri});
    ylim([0 1]);
    ylabel(yname);
    % ------------------------
    % axis ------------------------------------
    axs = get(ax1);
    ax1.TickLength = [0 0];
    ax1.FontSize = 10;
    % ----------------------------------------------------

    % save -------------------------------------------
    savename = sprintf('%s/%s',pSave,msr);
    printfig(savename,pSave,'w',w,'h',h,'closefig',1);
    % -------------------------------------------
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





%% save data ++++++++++++
cd(pSave); save('data.mat','MWTDB','MWTSet');
% ----------------------

















