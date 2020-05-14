function raster_tapResponseShift(start,finish,pMWT,varargin)
%% RL201511041222_dose_raster_stats_tapResponse(start, finish,pHome)
% Analyze after tap response from raster plot data
% Input
%     pHome = '/Users/connylin/Dropbox/Lab/MWT_Analysis_ByExp/Dose_10sISI';
%     start = start time by seconds i.e. 95; 
%     finish = end time by seconds, i.e. 105;

%% HYPOTHESIS: 
% reversal duration: all doses shortens reversal duration
% response latency: higher doses (400mM+) increases response latency
% reversal probability: probability started to decrease from 300mM < use
% bethoven
% separate data to:
%     1-after tap: + response latency
%     2-pre tap: use the same analysis as basal activity (not done)
%     after tap responses: changed within 0.4s (two frames) 
%     3-post tap response: same as basal activity




%% FUNCTION PATHS
% general home path
pMatFun = '/Users/connylin/Dropbox/Code/Matlab';
pRLFun = '/Users/connylin/Dropbox/RL/Code/Modules';
% add packges
addpath([pMatFun,'/General']);
addpath([pRLFun,'/Graphs']);
addpath([pRLFun,'/Graphs/rasterPlot_colorSpeed']);
%% DEFAULTS 
nInput = 3;
pAnalysis = '/Volumes/COBOLT/MWT';
analysisNLimit = 5;
alcoholTestDose = 400;
refStrain = 'N2';
tapNExpected = 30;
alphaValue = 0.05;
posttestname = 'bonferroni';
frameInt = 0.2;
NWORMS = Inf;
pSave = '/Users/connylin/Dropbox/RL/Dance Output';
%% VARARGIN PROCESSOR
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


%% PROCESS VARIABLES
pSave = sprintf('%s/%s',pSave,mfilename);
frameNumberExpected = numel((start:frameInt:finish-frameInt));

%% calculate raster plot with tap aligned in middle
fprintf('\n> Running raster plot analysis\n');
for x = 1:numel(pMWT)
    pMWTp = pMWT{x};
    trinitySummary(pMWTp,start,finish,'NWORMS',NWORMS,'frameInt',frameInt);
    return
end

%% tap response: summarize data for post tap
pSave = sprintf('%s/Stats/tapResponse_%d_%d',pHome,start,finish);
if isdir(pSave) == 0; mkdir(pSave); end
[~,~,~,pG] = dircontent(pHome);
pG(regexpcellout(pG,'(graffle|Stats)')) = []; % get rid of graffle files
rasterName = sprintf('*%d_%d*rasterData',start, finish);
DataSum = struct;
DataProb = struct;

for gi = 1:numel(pG)
    % get dir to MWT
    [~,fG] = fileparts(pG{gi});
    fprintf('Summarizing group %s\n',fG);
    i = find(ismember(plateInfo.groupname,fG));
    pMWT = cell(numel(i),1);
    A = plateInfo(i,:);
    for x = 1:numel(i)
        p = sprintf('%s/%s/%s/%s',pAnalysis,A.expname{x},A.groupname{x}, A.mwtname{x});
        pMWT{x} = p;
    end
    
   % load raster plot data
    B = table;
    D = []; FT = [];
    for mwti = 1:numel(pMWT)
        pMWTp = pMWT{mwti};
        fn = char(dircontent(pMWTp,sprintf('%ds_%ds_fint_%.1f_*_rasterdata.trinitySummary',start,finish,frameInt)));
        if isempty(regexp(fn,'(N0)','once')) == 1 % if N is larger than 1
            a = dlmread(fn);
            b = dlmread(char(dircontent(pMWTp,sprintf('%ds_%ds_fint_%.1f_*_frametime.trinitySummary',start,finish,frameInt))));
            c = dlmread(char(dircontent(pMWTp,sprintf('%ds_%ds_fint_%.1f_*_wormID.trinitySummary',start,finish,frameInt))));
            % if a is bigger than D, which happens when tap lands on a integer such as 100
            if size(a,2) > frameNumberExpected
                a = a(:,1:frameNumberExpected);
                b = b(:,1:frameNumberExpected);
            end
            if isempty(D) == 1 && isempty(a) == 0
                 D = a;
                 FT = b;
            elseif isempty(D) == 0
                 D = [D;a];
                 FT = [FT;b];
            end
            % create plate summary
            e = table;
            [~,fn] = fileparts(pMWT{mwti});
            e.fMWT = repmat(fn,size(a,1),1);
            e.wormid = c;
            B = [B;e];
        end
    end
    
    
    %% analyze after tap responses
    tapCol = (size(FT,2)/2)+1; % get tap column position
    % calculate onset and end of response shift before after tap
    tDirChange = nan(size(D,1),4);
    d = D(:,tapCol-1:end);
    d(d > 0) = 1;
    d(d < 0) = -1;
    d(d == 0) = 0;
    dDir = d;
    tDirChange(:,3) = dDir(:,1);
    for x = 2:size(d,2)
        d(d(:,1) ~= d(:,x),x) = nan;
    end
    for x = 1:size(d,1)
        a = find(isnan(d(x,:)));
        if isempty(a) == 0
            % find start time
            t1 = tapCol-2 + a(1);
            if isempty(t1) == 0
                tDirChange(x,1) = t1; 
                tDirChange(x,4) = dDir(x,a(1));
                
            end
            % find end time
            i = find(diff(a) > 1);
            if isempty(i) == 0 && numel(diff(a)) > 1
                tDirChange(x,2)  = tapCol-2 +(a(i(1))-1)+1;
            elseif isempty(i) == 1 && numel(a) == 1
                tDirChange(x,2)  = tapCol-2 +(a(1)-1)+1;
            end
        end
    end
    

    
    %% hypthesis 1: reversal duration: all doses shortens reversal duration
    fname = 'ResponseShiftDur';
    b = tDirChange(:,2); % end of frame
    b(isnan(b)) = frameNumberExpected;
    a = (b- tDirChange(:,1) +1).*frameInt;
    a(isnan(tDirChange(:,1))) = 0;
    dlmwrite(sprintf('%s/%s_%s.csv',pSave,fname,fG),a);
    DataSum.(fname).(fG) = a;
    
    %% hypothesis 3: response latency: higher doses (400mM+) increases response latency
    fname = 'ResponseShiftLatency';
    % calculate response latency (not very accurate becauese frames were summarized in 0.2 intervals
    a = (tDirChange(:,1) - tapCol+1).*frameInt;
    a(isnan(tDirChange(:,1))) = inf; % record no change as inf
    dlmwrite(sprintf('%s/%s_%s.csv',pSave,fname,fG),a);
    a(isnan(tDirChange(:,1))) = nan; % record no change as inf
    DataSum.(fname).(fG) = a;

    
    %% hypothesis 2: reversal probability: probability started to decrease from 300mM 
%     fname = 'ResponseRevProb';
%     % calculate response reversal probaiblity
%     d = D(:,tapCol-1:end);
%     % select worms responded within 1 second
%     k = tDirChange(:,1) <= ((frameNumberExpected/2)+1) + (1/frameInt);
%     % select worms that were pausing or forwarding before tap
%     i = tDirChange(:,3) ~= -1;
%     i = (k == true & i == true);
%     %%
%     G = B(i,1);
%     X = tDirChange(i,4) == -1;
%     mwtnames = unique(G);
%     for x = 1:numel
%         
%     end
%     %% select worms that reversed
%     j = 
%     %%
% 
%     %%
%     
%     DataProb.(fname) = table;
%     G = B(k,:);
%     a = 
%     X = 
%     sum(j)/sum(i)
%     sum(k)
% %     % calculate response latency (not very accurate becauese frames were summarized in 0.2 intervals
% %     a = (tDirChange(:,1) - tapCol+1).*frameInt;
% %     a(isnan(tDirChange(:,1))) = inf; % record no change as inf
% %     dlmwrite(sprintf('%s/%s_%s.csv',pSave,fname,fG),rDuration);
% %     a(isnan(tDirChange(:,1))) = nan; % record no change as inf
% %     DataSum.(fname).(fG) = a;
%     
    
end



%% stats
rNames = fieldnames(DataSum);
fid = fopen([pSave,'/anova.txt'],'w');
fprintf(fid,'ANOVA stats\n');
for ri = 1:numel(rNames)
    rNamesT = (rNames{ri});
    gnames = fieldnames(DataSum.(rNamesT));
    D = []; G = {};
    for g = 1:numel(gnames)
        d = DataSum.(rNamesT).(gnames{g});
        D = [D;d];
        G = [G; repmat(gnames(g),size(d,1),1)];  
    end
    % descriptive stats
    [gn,n,m,se] = grpstats(D,G,{'gname','numel','mean','sem'});
    T = table;
    T.groupname = gn;
    T.N = n;
    T.mean = m;
    T.SE = se;
    writetable(T,sprintf('%s/%s_stats.csv',pSave,rNamesT));
    % anova
    [p,t,s] = anova1(D,G,'off');
    fprintf(fid,'%s, %s\n',rNamesT,anova_textresult(t));
    [c,m2,~,~] = multcompare(s,'ctype',posttestname,'alpha',alphaValue, 'display', 'off');
    % determine pairs 
    d = c(:,3:5);
    d(d > 0) = 1;   
    d(d<0) = -1;
    i = c(abs(sum(d,2)) == 3,1:2);
    % export
    A = table;
    A.pair1 = gn(i(:,1));
    A.pair2 = gn(i(:,2));
    writetable(A,sprintf('%s/%s_%s_alpha_%0.1f.txt',pSave,rNamesT,posttestname,alphaValue*100));
end




%% report done
fprintf('* %s  done *\n',mfilename)










return








% 
% 
% 
% return
%     %%
%     
%     % average each worm's acceleration speed
%     varname = {'N_total','N_for','mean_for','SD_for','N_pause','N_rev',...
%         'mean_rev','SD_rev','N_revDurOnset','N_dirShift','mean_revDur','SD_revDur','mean_speed','SD_speed',...
%         'mean_speed_unbiased','SD_speed_unbiased'};
%     D = nan(size(a,1),numel(varname));
%     D = array2table(D,'VariableNames',varname);
% 
%     for wormi = 1:size(a,1)
%         D.N_total(wormi) = numel(a(wormi,:));
%         % forward speed
%         d = a(wormi,a(wormi,:) > 0);
%         D.N_for(wormi) = numel(d);
%         D.mean_for(wormi) = mean(d);
%         D.SD_for(wormi) = std(d);  
%         % get pause
%         D.N_pause(wormi) = numel(a(wormi,a(wormi,:) == 0));
%         % reversal
%         d = a(wormi,a(wormi,:) < 0);
%         D.N_rev(wormi) = numel(d);
%         D.mean_rev(wormi) = mean(d);
%         D.SD_rev(wormi) = std(d); 
%         % reversal duration
%         b = a(wormi,:);
%         i = find(b < 0)'; % get frame with reversal
%         if isempty(i) == 1
%             D.N_revDurOnset(wormi) = 0;
%         else
%             j = find([2;diff(i)]>1); % find first frame with reversal
%             D.N_revDurOnset(wormi) = numel(j);
%             dur = ([i([j(2:end)-1;numel(i)]) - i(j)] +1).*0.2; % calculate dur of continuous reversal
%             D.mean_revDur(wormi) = mean(dur);
%             D.SD_revDur(wormi) = std(dur);
%         end 
%         % speed 
%         D.mean_speed(wormi) = mean(a(wormi,:));
%         D.mean_speed(wormi) = std(a(wormi,:));
%         % speed (unbiased for dir, and no zeros)
%         d = a(wormi,:);
%         d(d == 0) = [];
%         d(d<0) = -d(d<0);
%         if isempty(d) == 1;
%             D.mean_speed_unbiased(wormi) = 0;
%             D.SD_speed_unbiased(wormi) = 0;
%         else
%             D.mean_speed_unbiased(wormi) = mean(d);
%             D.SD_speed_unbiased(wormi) = mean(d);
%         end
%         % direction shifting
%         d = a(wormi,:);
%         d(d>0) = 1;
%         d(d<0) = -1;
%         D.N_dirShift(wormi) = sum([0 diff(d)] ~= 0);
%         
%     end
%     writetable(D,sprintf('%s/%s_%d_%d.csv',pG{gi},prefixName,start, finish));
%     S{gi} = D;
%     
%     % calculate group mean
%     A.N_total(gi) = numel(D.mean_for); 
%     % forward speed
%     A.N_forward(gi) = sum(~isnan(D.mean_for)); % number of worms did any forward movement 
%     A.mean_forward(gi) = nanmean(D.mean_for);
%     A.SD_forward(gi) = nanstd(D.mean_for);
%     % forward freq
%     A.mean_forwardFreq(gi) = mean(D.N_for./D.N_total);
%     A.SD_forwardFreq(gi) = std(D.N_for./D.N_total);
%     % proportion paused 
%     d = D.N_pause/mean(D.N_total);
%     A.mean_pauseFreq(gi) = mean(d);
%     A.SD_pauseFreq(gi) = std(d);
%     % proportion rev
%     d = D.N_rev./D.N_total;
%     A.mean_revFreq(gi) = mean(d);
%     A.SD_revFreq(gi) = std(d);
%     % reversal speed
%     A.N_rev(gi) = sum(~isnan(D.mean_rev));
%     A.mean_rev(gi) = nanmean(D.mean_rev);
%     A.SD_rev(gi) = nanstd(D.mean_rev);
%     % reversal duration
%     A.mean_revDur(gi) = nanmean(D.mean_revDur);
%     A.SD_revDur(gi) = nanstd(D.mean_revDur);
%     % speed 
%     A.mean_speed(gi) = mean(D.mean_speed);
%     A.SD_speed(gi) = std(D.mean_speed);
%     % speed (unbiased for dir, and no zeros)
%     A.mean_speed_unbiased(gi)  = mean(D.mean_speed_unbiased);
%     A.SD_speed_unbiased(gi) = std(D.SD_speed_unbiased);
%     % direction shifting
%     A.mean_dirShift(gi) = mean(D.N_dirShift);
%     A.SD_dirShift(gi) = std(D.N_dirShift);
% end
% 
% % create output table;
% [~,gname] = cellfun(@fileparts,pG,'UniformOutput',0);
% T = table;
% T.gname = gname;
% T = [T A];
% writetable(T,sprintf('%s/%s_%s_%s_summary.csv',pSave,prefixName,start, finish));
% 
% 
% 
% % statistics
% pSave = sprintf('%s/Stats/%s_%d_%d_stats',pHome,prefixName,start, finish);
% if isdir(pSave) == 0; mkdir(pSave); end
% varname = {'N_total',...
% 'speed','forward_speed','forward_freq','pause_freq',...
% 'rev_freq','rev_speed','rev_dur','speed_unbiased',...
% 'dirShift'};
% 
% alphaValue = 0.05;
% posttestname = 'bonferroni';
% fid = fopen([pSave,'/anova.txt'],'w');
% fprintf(fid,'ANOVA stats\n');
% for v = 1:numel(varname)
%     
%     vT = varname{v};
%     G = []; D = [];
%     
%     switch vT
%         case 'N_total'
%             for g =1:numel(pG)
%                 d = S{g}.mean_for;
%                 D = [D;d]; G = [G; repmat(g,numel(d),1)];
%             end
%         case 'forward_freq' % percent of worms moved forward / total worms
%             for g =1:numel(pG)
%                 d = S{g}.N_for./S{g}.N_total;
%                 D = [D;d]; G = [G; repmat(g,numel(d),1)];
%             end
%         case 'rev_freq' % perecen tof worms reversed/total worms
%             for g =1:numel(pG)
%                 d = S{g}.N_rev./S{g}.N_total;
%                 D = [D;d]; G = [G; repmat(g,numel(d),1)];
%             end
%         case 'speed' 
%             for g =1:numel(pG)
%                 d = S{g}.mean_speed;
%                 D = [D;d]; G = [G; repmat(g,numel(d),1)];
%             end
%         case 'forward_speed' 
%             for g =1:numel(pG)
%                 d = S{g}.mean_for;
%                 D = [D;d]; G = [G; repmat(g,numel(d),1)];
%             end
%         case 'pause_freq' 
%             for g =1:numel(pG)
%                 d = S{g}.N_pause./S{g}.N_total;
%                 D = [D;d]; G = [G; repmat(g,numel(d),1)];
%             end
%         case 'rev_speed' 
%             for g =1:numel(pG)
%                 d = S{g}.mean_rev;
%                 D = [D;d]; G = [G; repmat(g,numel(d),1)];
%             end
%         case 'rev_dur' 
%             for g =1:numel(pG)
%                 d = S{g}.mean_revDur;
%                 D = [D;d]; G = [G; repmat(g,numel(d),1)];
%             end
%         case 'speed_unbiased' 
%             for g =1:numel(pG)
%                 d = S{g}.mean_speed_unbiased;
%                 D = [D;d]; G = [G; repmat(g,numel(d),1)];
%             end
%         case 'dirShift' 
%             for g =1:numel(pG)
%                 d = S{g}.N_dirShift;
%                 D = [D;d]; G = [G; repmat(g,numel(d),1)];
%             end
%     end
%     % anova
%     [p,t,s] = anova1(D,G,'off');
%     fprintf(fid,'%s, %s\n',vT,anova_textresult(t));
%     [c,m,~,~] = multcompare(s,'ctype',posttestname,'alpha',alphaValue, 'display', 'off');
%     % determine pairs 
%     d = c(:,3:5);
%     d(d > 0) = 1;
%     d(d<0) = -1;
%     i = c(abs(sum(d,2)) == 3,1:2);
%     % export
%     A = table;
%     A.pair1 = gname(i(:,1));
%     A.pair2 = gname(i(:,2));
%     writetable(A,sprintf('%s/%s_%s_alpha_%0.1f.txt',pSave,vT,posttestname,alphaValue*100));
%     
% end
% fclose(fid);














