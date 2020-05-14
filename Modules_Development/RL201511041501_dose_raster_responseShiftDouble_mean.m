%% HYPOTHESIS: 
%     after tap, (at later tap) worms do forward > short reverse > pause
%     alcohol knocks out this resposne
%     find, for worms who does pause/forward > reverse > pause 
%     pause/forward > reverse (mean duration) > pause (mean duration)



%% paths & settings
pHome = '/Users/connylin/Dropbox/Lab/MWT_Analysis_ByExp/Dose_10sISI';
start = 385; 
finish = 395;
taptime = 390;
frameInt = 0.1;
NWORMS = Inf;
tapCol = ((taptime-start)/frameInt)+1;
alphaValue = 0.05;
posttestname = 'bonferroni';
frameNumberExpected = numel((start:frameInt:finish-frameInt));

% exp info path
pSave = sprintf('%s/Stats/tapResponse_%d_%d',pHome,start,finish);
if isdir(pSave) == 0; mkdir(pSave); end
pAnalysis = '/Users/connylin/Dropbox/Lab/MWT_Analysis';
% function path
addpath('/Users/connylin/Dropbox/MATLAB/Function_Library_Public');
addpath('/Users/connylin/Dropbox/MATLAB/programs_rankinlab/Library/Modules/rasterPlot_colorSpeed');


%% calculation
% load global data
plateInfo = readtable(sprintf('%s/plate_info.csv',pHome));
[~,~,~,pG] = dircontent(pHome);
pG(regexpcellout(pG,'(graffle|Stats)')) = []; % get rid of graffle files
Data = {};
gname = {};
for gi = 1:numel(pG)
    %% get dir to MWT
    [~,fG] = fileparts(pG{gi});
    gname{gi}= fG;
    fprintf('Summarizing group %s\n',fG);
    i = find(ismember(plateInfo.groupname,fG));
    pMWT = cell(numel(i),1);
    A = plateInfo(i,:);
    for x = 1:numel(i)
        p = sprintf('%s/%s/%s/%s',pAnalysis,A.expname{x},A.groupname{x}, A.mwtname{x});
        pMWT{x} = p;
    end
    
    %% get trinity data
    for mwti = 1:numel(pMWT)
        pMWTp = pMWT{mwti};
        addpath('/Users/connylin/Dropbox/MATLAB/programs_rankinlab/Library/Modules/rasterPlot_colorSpeed');
        fprintf('\n> Running raster plot analysis\n');
        a = trinitySummary(pMWTp,start,finish,'frameInt',frameInt,'NWORM',Inf,'saveoption',0);        
        % if a is bigger than D, which happens when tap lands on a integer such as 100
        if size(a,2) > frameNumberExpected
            a = a(:,1:frameNumberExpected);
        end
        if mwti == 1; D = a; else D = [D;a]; end
    end
    Data{gi} = D;
    %% mean response rate
    if gi == 1
        M = mean(D);
        N = size(D,1);
        SE = std(D)./sqrt(size(D,1)-1);
    else
        M = [M;mean(D)];
        N = [N;size(D,1)];
        SE = [SE;std(D)./sqrt(size(D,1)-1)];
    end
end


%% plot
T = table;
T.time = [start:frameInt:finish-frameInt]';
a = array2table(M','VariableNames',gname);
T = [T a];
v = cellfun(@strcat,gname',cellfunexpr(gname','_SE'),'UniformOutput',0);
b = array2table(SE','VariableNames',v);
T = [T b];
fn = sprintf('%s/MeanSpeed_descriptive_stats.csv',pSave);
writetable(T,fn);


return
    