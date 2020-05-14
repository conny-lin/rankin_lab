%% INPUT
pData = '/Volumes/COBOLT/MWT';
pDest = '/Users/connylin/Dropbox/RL/MWT_Analysis_ByExp/STH/10sISI/by strain rasterNew';
% create time frames
startTime = 88;
endTime = 405;
int = 10;
assaywindow = 30;
frameInt = 0.2;
NWORMS = Inf;
expLimit = 'N2 within';
%% create time frame
% startList = startTime:int:endTime;
startList = [98 238 368];
endList = startList+assaywindow;
% nAssayFrame = assaywindow/frameInt;

%% FUNCTION PATHS
% general home path
pMatFun = '/Users/connylin/Dropbox/Code/Matlab';
pRLFun = '/Users/connylin/Dropbox/RL/Code/Modules';
% add packges
addpath([pRLFun,'/Graphs/rasterPlot_colorSpeed']);
addpath([pRLFun,'/MWTDatabase']);
addpath([pRLFun,'/Chor']);
addpath([pMatFun,'/General']);
addpath([pMatFun,'/Graphs']);


%% construct strain list
p = '/Users/connylin/Dropbox/RL/MWT_Analysis_ByExp/STH/10sISI/by strain vs N2/Strains';
[~,~,strainList] = dircontent(p);

%% run
% strainSuccess = true;
strainList = {'CB1112','KG571','JPS428','KP4','VM3109','JPS344','JPS338','JPS326','NM1815','RB781',};
for si = 1:numel(strainList)
    strainName = strainList{si};
% strainName = 'VG301';

% construct group names
% groupnameList = [{'N2','N2_400mM'} strainName {[char(strainName),'_400mM']}];
    searchgroup = [strainName {[char(strainName),'_400mM']}];
    ctrlgroup = {'N2','N2_400mM'};
    groupnameList = [searchgroup ctrlgroup];
% groupnameList = {'N2','N2_400mM'};


% query data
%     DbT = MWTDatabase_query(pData,'groupname',groupnameList,...
%             'gnameSearchType','within','rc','100s30x10s10s');
% LOAD DATABASE
D = load([pData,'/MWTDatabase.mat']);
MWTDatabase = D.MWTDatabase;
% get targets
Db = MWTDatabase.mwt;
i = ismember(Db.groupname,searchgroup) &...
    ismember(Db.rc,'100s30x10s10s');
i = ismember(Db.expname,unique(Db.expname(i))) & ...
    ismember(Db.groupname,groupnameList);
DbT = Db(i,:);

% exclude bad exp (default)
badExp = {'20140906C_NG_100s30x10s10s_npr','20130611C_NG_100s30x10s10s_100Scannotfix'};
DbT(ismember(DbT.expname,badExp),:) = [];

unique(DbT.expname)
unique(DbT.groupname)
size(DbT)


% run raster plot
for gi = 1:numel(groupnameList)
    gT = groupnameList{gi};
    pMWT = DbT.mwtpath(ismember(DbT.groupname,gT));
    if isempty(pMWT)== 0

        % chor
        pMWTc = convertTrinityDat2Mat(pMWT,1); 
        L = chormaster4('Trinity',pMWTc);
        % summarize trinity data and delete .dat file to save memory
        pMWTbad = convertTrinityDat2Mat(pMWTc,1); 
        % exclude bad files
        pMWToriginal = pMWT;
        pMWT(ismember(pMWT,pMWTbad)) = [];

        % raster plot per time points
        for ti = 1:numel(startList)
            % create save folder
            pSave = sprintf('%s/%s/%s/%s',pDest,strainName,gT);
            if isdir(pSave) == 0; mkdir(pSave); end

            % time
    %                 start = 268;
    %                 finish = 308;
            start = startList(ti);
            finish = endList(ti);
    %                 
            % run raster
            [f1,Data,savename,Tog,mwtid,rTime,Import] = rasterPlot_colorSpeed(pMWT,start,finish,...
                'NWORMS',NWORMS,'visibleG',0,'frameInt',frameInt);
            % save trinity import
            if isempty(dircontent(pSave,'trinity.mat')) == 1
                cd(pSave); save('trinity.mat','Import');
            end
            % save fig
            cd(pSave);
            set(f1,'PaperPositionMode','auto'); % set to save as appeared on screen
            print (f1,'-depsc', '-r1200', savename); % save as eps
            close;
            % save data
            save(sprintf('%s/%s rasterData.mat',pSave,savename),'pMWT','Data','mwtid','Tog','rTime');
        end
    end
end


%% electrophys graph
for ti = 1:numel(startList)
t1 = startList(ti);
t2 = endList(ti);
T = table;
for gi = 1:numel(groupnameList)
    gT = groupnameList{gi};
    pSave = sprintf('%s/%s/%s/%s',pDest,strainName,gT);
    pSave = regexprep(pSave,'[/]\>','');
    if isdir(pSave) == 1
    [f,p] = dircontent(pSave);
    if isempty(f) == 0
        str = sprintf('[_]%d[_]%d(?=[_])',t1,t2);
        i = regexpcellout(f,str) & regexpcellout(f,'mat\>');

        D = load(p{i},'Data');
        D = D.Data;
        M = mean(D);
    %     N = size(D,1);
    %     SE = std(D)./sqrt(size(D,1)-1);
        T.(gT) = M';
    end
    end
end

pHome =fileparts(pSave);
writetable(T,sprintf('%s/ef graph %d-%d.csv',pHome,t1,t2));
end

end

fprintf('\n\nDONE, and below strains did not work:\n');
% disp(char(strainList(~strainSuccess)));






return
% D = Data;
% % plot all data
% plot(D','Color',[0.8 0.8 .8]); hold on
% savefigepsOnly150([savename,' ef all'],pSave);
% 
% M = mean(D);
% N = size(D,1);
% SE = std(D)./sqrt(size(D,1)-1);
% 
% % plot(M+SE,'Color',[0 0 0])
% t = table; 
% t.mean = M';
% t.mse1 = (M+SE)';
% t.mse2 = (M-SE)';
% writetable(t,sprintf('%s/ef graph.csv',pSave));
% 
% %%
% 
% 
% bv = min(M-(2*SE));
% area(M+(2*SE),'basevalue',bv,'FaceColor',[.5 .5 .5],'EdgeColor','none')
% hold on
% 
% area(M-(2*SE),'basevalue',bv,'FaceColor',[1 1 1],'EdgeColor','none')
% % plot(M-SE,'Color',[0 0 0])
% line([1:size(D,2)+1],zeros(1,size(D,2)+1),'Color',[0 0 0])




%% TROEBLE SHOOTING - remove no tap plate
% load raster plot .mat output
cd(pSave); 
load('rasterPlot_98_108_N_680 rasterData.mat','Data');
load('rasterPlot_98_108_N_680 rasterData.mat','rTime');
load('rasterPlot_98_108_N_680 rasterData.mat','mwtid');
%% find no reversals
t1 = 11;
t2 = 15;
d = Data(:,t1:t2);
rev = sum(d < 0.2,2); % reversals
pause = sum(d==0,2); % pause
[s,gn,n] = grpstats(rev,mwtid,{'sum','gname','numel'});
t = table;
t.mwtid = cellfun(@str2num,gn);
t.sum = s;
t.N = n.*(t2-t1+1);
t.p = t.sum./t.N;
t
% selection
ip = t.mwtid(t.p > 0.8);



A = parseMWTinfo(pMWT);
A.expname(ip)




























