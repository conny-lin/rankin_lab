function meanspeed(mutant,pSaveHome)


%% addpath
addpath('/Users/connylin/Dropbox/rl/Code/Modules/Graphs/rasterPlot_colorSpeed');
addpath('/Users/connylin/Dropbox/Code/Matlab/Library/General');
addpath_allsubfolders('/Users/connylin/Dropbox/Code/Matlab/Library/General');
addpath('/Users/connylin/Dropbox/rl/Code/Modules/Chor');
addpath('/Users/connylin/Dropbox/RL/Code/Modules/MWTDatabase');
addpath('/Users/connylin/Dropbox/RL/Code/Modules/Stats');
pData = '/Volumes/COBOLT/MWT';

%% setting
% time setting
runcond = '100s30x10s10s';
expectedTaptime = [100:10:(100+(10*29))];
assaystarttime = 10;
assayendtime = expectedTaptime(end)+10;
speed_unit = 0.1;
% condition setting
wildtype = 'N2';
cond = {'400mM'};
% 
% 
% 
%% process setting variables
% time
expectedTapN = numel(expectedTaptime);
% make groupname list
strain = [{wildtype} {mutant}];
nrow = numel(strain)*(numel(cond)+1);
groupnameList = cell(nrow,1);
n = 1;
for x = 1:numel(strain)
   groupnameList(n) = strain(x);
   n = n+1;
   groupnameList{n} = strjoin([strain(x),{'_'},cond],'');
   n = n+1;
end
% 
% 
% 
%% get pMWT
% load database
load('/Users/connylin/Dropbox/rl/MWTDB/MWTDB.mat')
Database = MWTDB.text; clear MWTDB;
% find experiment names
i = ismember(Database.strain, mutant) ...
    & ismember(Database.rc, runcond) ...
    & Database.exp_date > 20111213; 
j = ismember(Database.expname, unique(Database.expname(i)));
DBT = Database(j,:);
DBT(~ismember(DBT.groupname, groupnameList),:) = [];
fprintf('\ngroupname:\n');
disp(char(unique(DBT.groupname)))
fprintf('\nfrom exp:\n');
disp(char(unique(DBT.expname)))
pMWT = DBT.mwtpath;
clear Database;
fprintf('\nDatabase found %d MWT folders\n',numel(pMWT));
% 
% 
% 
%% CHOR
pMWTc = convertTrinityDat2Mat(pMWT,1); 
L = chormaster4('Trinity',pMWTc);
% summarize trinity data and delete .dat file to save memory
pMWTbad = convertTrinityDat2Mat(pMWTc,1); 
% exclude bad files
pMWToriginal = pMWT;
pMWT(ismember(pMWT,pMWTbad)) = [];
Db = parseMWTinfo(pMWT);
% add problem output col
Db.problem = cell(size(Db,1),1);
groupnameList = unique(Db.groupname);
fprintf('\n%d MWT folders contain valid chor files\n',size(Db,1));
% 
% 
% 
%% DEALING WITH TRINITY DATA ==============================================
% load trinity data
Trinity = cell(size(pMWT));
for mwti = 1:numel(pMWT)
    fprintf('loading %d/%d MWT plate\n',mwti,numel(pMWT));
    % get trinity.mat
    pMWTp = pMWT{mwti};
    Trinity{mwti} = load_Trinity(pMWTp);
end
% process trinity data
SumDatap = cell(size(pMWT));
for mwti = 1:numel(pMWT);
    % find out total rows
    [nrow, ~] = cellfun(@size,Trinity{mwti}(:,2));
    nrowsum = sum(nrow);
    SumData = nan(nrowsum,5);
    rowstart = [0;cumsum(nrow)];
    % get only time(1), tap(7), speed (4), bias (6)
    wormid = Trinity{mwti}(:,1);
    wormid = cellfun(@str2num,wormid);
    for wrmi = 1:numel(wormid)
        wormid_current = wormid(wrmi);
        Data = Trinity{mwti}{wrmi,2};
        speed = Data(:,4);
        bias = Data(:,6);
        speed_dir = speed.*bias;
        time = Data(:,1);
        tap = Data(:,7);
        wid = repmat(wormid_current,size(Data,1),1);
        pid = repmat(mwti,size(Data,1),1);        
        D = [pid wid time tap speed_dir];
        SumData(rowstart(wrmi)+1:rowstart(wrmi+1),1:5) = D;
    end
    SumDatap{mwti} = SumData;
end
% put data in D
D = cell2mat(SumDatap);
clear Trinity SumDatap SumData;
% convert to table
D = array2table(D,'VariableNames',{'mwtid','wormid','time','tap','speed_bias'});
clear RawData;
% 
% 
% 
%% VALIDATIONS  ========================================================
fprintf('keep time between %d-%ds\n',assaystarttime,assayendtime);
% remove time before assay time
D(D.time < assaystarttime,:) = [];
% remove time after assay time
D(D.time > assayendtime,:) = [];

% adjust time according to tap time
itap = D.tap > 0; % index to tap time
Taptime = [D.mwtid(itap) D.time(itap) round(D.time(itap))]; % get tap time info
% remove plates with wrong tap times
taptime_round = Taptime(:,3);
i = ~ismember(taptime_round,expectedTaptime);
mwtid = Taptime(i,1);
if isempty(i) == 0
    mwtid_wrongtap = unique(mwtid);
    % report info from database
    j = ismember(Db.mwtid, mwtid_wrongtap);
    fprintf('\nremove plates below has incorrect tap time:\n');
    disp([Db.expname(j) Db.mwtname(j) Db.groupname(j)]);
    Db.problem(j) = {'incorrect tap time'};
    % delete from raw data
    D(ismember(D.mwtid, mwtid_wrongtap),:) = [];
end


%% create legend for tap time
itap = D.tap > 0; % index to tap time
Taptime = [D.mwtid(itap) D.time(itap) round(D.time(itap))]; % get tap time info
tpu = unique(Taptime(:,3)); % find unique tap time
ntap = numel(tpu); % tap time number
fprintf('\n%d unique tap times\n',ntap);
for n = 1:numel(tpu)
    i = Taptime(:,3) == tpu(n);
    Taptime(i,3) = n; 
end
% find tap time range per plate
Tapstart = nan(ntap,numel(pMWT));
Tapend = Tapstart;
for n = 1:numel(tpu)
    i = Taptime(:,3) == n;
    mwt = Taptime(i,1);
    time =Taptime(i,2);
    [m1,m2,mwtid] = grpstats(time,mwt,{'min','max','gname'});
    plateid = cellfun(@str2num,mwtid);
    Tapstart(n,plateid) = m1';
    Tapend(n,plateid) = m2';
end



%% exclude plates with missing taps
badplateid = any(isnan(Tapstart)) | any(isnan(Tapend));
if sum(badplateid) > 1
    fprintf('plates below do not have all taps:\n');
    disp([Db.expname(badplateid) Db.groupname(badplateid) Db.mwtname(badplateid)]);
    fprintf('remove from analysis\n');
    D(ismember(D.mwtid,find(badplateid)),:) = [];
end


%% standardize time
% use the first plate as standard
taptimeStd = expectedTaptime(1);
dt = Tapstart(1,:) - taptimeStd;
% adjust time
D.time_adj = D.time + dt(D.mwtid)';
% create time by 0.02 seconds
D.time_adj_sec = round(D.time_adj./speed_unit).*speed_unit;

% save raw data
RawData = D;
cd(pSaveHome); save('data.mat','RawData','-v7.3');

%% create group id
groupname = unique(Db.groupname);
groupnameInd = table;
groupnameInd.group_id = (1:numel(groupname))';
groupnameInd.groupname = groupname;
[~,i] = ismember(Db.groupname,groupname);
D.group_id = i(D.mwtid);


%% calculate speed 2
fprintf('\ngroup speed data by group\n');
% row = worm, col = time
% convert to matrix
tu = unique(D.time_adj_sec);
ByGroup = struct;
ByGroup.time = tu;
for gi = 1:numel(groupname)
    fprintf('group %d/%d: %s\n',gi,numel(groupname),groupname{gi});
    i = D.group_id == gi;
    t = D.time_adj_sec(i);
    s = D.speed_bias(i);
    % convert raw data
    a = tabulate(t);
    B = nan(max(a(:,2)),numel(tu)); % output array
    for ti = 1:numel(tu)
        j = t == tu(ti);
        d = s(j);
        n = numel(d);
        B(1:n,ti) = d;
    end
    ByGroup.Speed.(groupname{gi}) = B;
end

%% calculate stats
fprintf('\ncalculate speed stats\n');
Stats = struct;
tu = ByGroup.time;
D = ByGroup.Speed;
for gi = 1:numel(groupname)
    B = D.(groupname{gi});
    A =  stats_descriptive_matrix(B,tu,'time');
    Stats.(groupname{gi}) = A;
end



%% export
fprintf('\nexporting matlab files\n');
cd(pSaveHome); 
save('data.mat','ByGroup','Stats','Db','Tapstart','Tapend','wildtype','mutant','cond','groupnameList','-append');



%% export to csv
fprintf('\nexporting csv files\n');
cd(pSaveHome);
D = Stats;
for gi = 1:numel(groupname)
    gn = groupname{gi};
    B = D.(gn);
    writetable(B,[gn,'.csv']);
end

%% export summary to csv
fprintf('\nexporting csv files\n');
cd(pSaveHome);
D = Stats;
nrow = numel(ByGroup.time);
ncol = numel(groupnameList);
B = nan(nrow,ncol);
for gi = 1:numel(groupnameList)
    gn = groupnameList{gi};
    B(:,gi) = D.(gn).mean;
end
tu = ByGroup.time;
B = array2table([tu B],'VariableNames',[{'time'};groupnameList]);
writetable(B,'mean_speed.csv');



%% plot - mean speed
fprintf('\nmaking graphs\n');
D = Stats;
t1 = D.(wildtype).time;
y1 = D.(wildtype).mean;
t2 = D.([wildtype,'_400mM']).time;
y2 = D.([wildtype,'_400mM']).mean;
figure1 = figure;
plot(t1,y1,'Color','k'); hold on
plot(t2,y2,'Color','r');
cd(pSaveHome);
savefig(figure1,wildtype);
set(figure1,'PaperSize',[15 8.5])
print(figure1,wildtype,'-depsc');
close;

t1 = D.(mutant).time;
y1 = D.(mutant).mean;
t2 = D.([mutant,'_400mM']).time;
y2 = D.([mutant,'_400mM']).mean;
figure1 = figure;
plot(t1,y1,'Color','b'); hold on
plot(t2,y2,'Color','m');
cd(pSaveHome);
savefig(figure1,mutant);
print(figure1,mutant,'-depsc');
close;

%% report done 
fprintf('\nDONE\n');




























































