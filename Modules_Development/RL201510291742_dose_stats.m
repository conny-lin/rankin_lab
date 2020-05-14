%% paths
% exp info path
pHome = '/Users/connylin/Dropbox/Lab/MWT_Analysis_ByExp/Dose_10sISI';
pSave = [pHome,'/Stats'];
% function path
addpath('/Users/connylin/Dropbox/MATLAB/Function_Library_Public');
pAnalysis = '/Users/connylin/Dropbox/Lab/MWT_Analysis';


%% HYPOTHESIS 1: 
% initial responses - faster acceleration speed in lower doses (100-300mM),
% slower acceleration speed in higher doses (400mM+)
% higher spontaneous reversal

%% get data
[~,~,~,pG] = dircontent(pHome);
pG(regexpcellout(pG,'(graffle|Stats)')) = []; % get rid of graffle files

varnameA = {'N_total','N_forward','N_rev',...
        'mean_speed',...
        'mean_forward','mean_forwardFreq','mean_pauseFreq','mean_revFreq',...
        'mean_rev','mean_revDur','mean_speed_unbiased','mean_dirShift',...
        'SD_speed','SD_forward','SD_forwardFreq','SD_pauseFreq','SD_revFreq','SD_rev',...
        'SD_revDur','SD_speed_unbiased','SD_dirShift'};
A = nan(size(pG,1),numel(varnameA));
A = array2table(A,'VariableNames',varnameA);
S = cell(size(pG));
for gi = 1:numel(pG)
    % get pre plate data
    [~,pRaster] = dircontent(pG{gi},'*85_95*rasterData');
    a = dlmread(pRaster{ri});
    a(1,:) = []; % eliminate first row
    
    % average each worm's acceleration speed
    varname = {'N_total','N_for','mean_for','SD_for','N_pause','N_rev',...
        'mean_rev','SD_rev','N_revDurOnset','N_dirShift','mean_revDur','SD_revDur','mean_speed','SD_speed',...
        'mean_speed_unbiased','SD_speed_unbiased'};
    D = nan(size(a,1),numel(varname));
    D = array2table(D,'VariableNames',varname);

    for wormi = 1:size(a,1)
        D.N_total(wormi) = numel(a(wormi,:));
        % forward speed
        d = a(wormi,a(wormi,:) > 0);
        D.N_for(wormi) = numel(d);
        D.mean_for(wormi) = mean(d);
        D.SD_for(wormi) = std(d);  
        % get pause
        D.N_pause(wormi) = numel(a(wormi,a(wormi,:) == 0));
        % reversal
        d = a(wormi,a(wormi,:) < 0);
        D.N_rev(wormi) = numel(d);
        D.mean_rev(wormi) = mean(d);
        D.SD_rev(wormi) = std(d); 
        % reversal duration
        b = a(wormi,:);
        i = find(b < 0)'; % get frame with reversal
        if isempty(i) == 1
            D.N_revDurOnset(wormi) = 0;
        else
            j = find([2;diff(i)]>1); % find first frame with reversal
            D.N_revDurOnset(wormi) = numel(j);
            dur = ([i([j(2:end)-1;numel(i)]) - i(j)] +1).*0.2; % calculate dur of continuous reversal
            D.mean_revDur(wormi) = mean(dur);
            D.SD_revDur(wormi) = std(dur);
        end 
        % speed 
        D.mean_speed(wormi) = mean(a(wormi,:));
        D.mean_speed(wormi) = std(a(wormi,:));
        % speed (unbiased for dir, and no zeros)
        d = a(wormi,:);
        d(d == 0) = [];
        d(d<0) = -d(d<0);
        if isempty(d) == 1;
            D.mean_speed_unbiased(wormi) = 0;
            D.SD_speed_unbiased(wormi) = 0;
        else
            D.mean_speed_unbiased(wormi) = mean(d);
            D.SD_speed_unbiased(wormi) = mean(d);
        end
        % direction shifting
        d = a(wormi,:);
        d(d>0) = 1;
        d(d<0) = -1;
        D.N_dirShift(wormi) = sum([0 diff(d)] ~= 0);
        
    end
    writetable(D,sprintf('%s/speed_85_95.csv',pG{gi}));
    S{gi} = D;
    
    % calculate group mean
    A.N_total(gi) = numel(D.mean_for); 
    % forward speed
    A.N_forward(gi) = sum(~isnan(D.mean_for)); % number of worms did any forward movement 
    A.mean_forward(gi) = nanmean(D.mean_for);
    A.SD_forward(gi) = nanstd(D.mean_for);
    % forward freq
    A.mean_forwardFreq(gi) = mean(D.N_for./D.N_total);
    A.SD_forwardFreq(gi) = std(D.N_for./D.N_total);
    % proportion paused 
    d = D.N_pause/mean(D.N_total);
    A.mean_pauseFreq(gi) = mean(d);
    A.SD_pauseFreq(gi) = std(d);
    % proportion rev
    d = D.N_rev./D.N_total;
    A.mean_revFreq(gi) = mean(d);
    A.SD_revFreq(gi) = std(d);
    % reversal speed
    A.N_rev(gi) = sum(~isnan(D.mean_rev));
    A.mean_rev(gi) = nanmean(D.mean_rev);
    A.SD_rev(gi) = nanstd(D.mean_rev);
    % reversal duration
    A.mean_revDur(gi) = nanmean(D.mean_revDur);
    A.SD_revDur(gi) = nanstd(D.mean_revDur);
    % speed 
    A.mean_speed(gi) = mean(D.mean_speed);
    A.SD_speed(gi) = std(D.mean_speed);
    % speed (unbiased for dir, and no zeros)
    A.mean_speed_unbiased(gi)  = mean(D.mean_speed_unbiased);
    A.SD_speed_unbiased(gi) = std(D.SD_speed_unbiased);
    % direction shifting
    A.mean_dirShift(gi) = mean(D.N_dirShift);
    A.SD_dirShift(gi) = std(D.N_dirShift);
end

% create output table;
[~,gname] = cellfun(@fileparts,pG,'UniformOutput',0);
T = table;
T.gname = gname;
T = [T A];
writetable(T,sprintf('%s/speed_85_95_summary.csv',pSave));



%% statistics
pSave = [pHome,'/Stats/speed_85_95_stats'];
if isdir(pSave) == 0; mkdir(pSave); end
varname = {'N_total',...
'speed','forward_speed','forward_freq','pause_freq',...
'rev_freq','rev_speed','rev_dur','speed_unbiased',...
'dirShift'};

alphaValue = 0.05;
posttestname = 'bonferroni';
fid = fopen([pSave,'/anova.txt'],'w');
fprintf(fid,'ANOVA stats\n');
for v = 1:numel(varname)
    
    vT = varname{v};
    G = []; D = [];
    
    switch vT
        case 'N_total'
            for g =1:numel(pG)
                d = S{g}.mean_for;
                D = [D;d]; G = [G; repmat(g,numel(d),1)];
            end
        case 'forward_freq' % percent of worms moved forward / total worms
            for g =1:numel(pG)
                d = S{g}.N_for./S{g}.N_total;
                D = [D;d]; G = [G; repmat(g,numel(d),1)];
            end
        case 'rev_freq' % perecen tof worms reversed/total worms
            for g =1:numel(pG)
                d = S{g}.N_rev./S{g}.N_total;
                D = [D;d]; G = [G; repmat(g,numel(d),1)];
            end
        case 'speed' 
            for g =1:numel(pG)
                d = S{g}.mean_speed;
                D = [D;d]; G = [G; repmat(g,numel(d),1)];
            end
        case 'forward_speed' 
            for g =1:numel(pG)
                d = S{g}.mean_for;
                D = [D;d]; G = [G; repmat(g,numel(d),1)];
            end
        case 'pause_freq' 
            for g =1:numel(pG)
                d = S{g}.N_pause./S{g}.N_total;
                D = [D;d]; G = [G; repmat(g,numel(d),1)];
            end
        case 'rev_speed' 
            for g =1:numel(pG)
                d = S{g}.mean_rev;
                D = [D;d]; G = [G; repmat(g,numel(d),1)];
            end
        case 'rev_dur' 
            for g =1:numel(pG)
                d = S{g}.mean_revDur;
                D = [D;d]; G = [G; repmat(g,numel(d),1)];
            end
        case 'speed_unbiased' 
            for g =1:numel(pG)
                d = S{g}.mean_speed_unbiased;
                D = [D;d]; G = [G; repmat(g,numel(d),1)];
            end
        case 'dirShift' 
            for g =1:numel(pG)
                d = S{g}.N_dirShift;
                D = [D;d]; G = [G; repmat(g,numel(d),1)];
            end
    end
    % anova
    [p,t,s] = anova1(D,G,'off');
    fprintf(fid,'%s, %s\n',vT,anova_textresult(t));
    [c,m,~,~] = multcompare(s,'ctype',posttestname,'alpha',alphaValue, 'display', 'off');
    % determine pairs 
    d = c(:,3:5);
    d(d > 0) = 1;
    d(d<0) = -1;
    i = c(abs(sum(d,2)) == 3,1:2);
    % export
    A = table;
    A.pair1 = gname(i(:,1));
    A.pair2 = gname(i(:,2));
    writetable(A,sprintf('%s/%s_%s_alpha_%0.1f.txt',pSave,vT,posttestname,alphaValue*100));
    
end
fclose(fid);

return

    






%% report done
fprintf('* %s  done *\n',mfilename)



return



















