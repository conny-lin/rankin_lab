function [Stats, DataG] = ephys_stats(DataG,pSave)

%% ANALYSIS
Stats = struct;

%% Process all sample --------

% find n(rows) and time (columns)
samplesize = nan(size(DataG,2),1);
timepts = samplesize;
for gi =1:size(DataG,2)
    [DataG(gi).N, DataG(gi).timepts] = size(DataG(gi).speedbm);
end


% decide group sequence
gu = cell(size(DataG,2),1);
for gi =1:size(DataG,2)
    gu{gi} = DataG(gi).name;
end
guseq = [find(regexpcellout(gu,'N2'));find(~regexpcellout(gu,'N2'))]';

% graph mean curve 
color = [0 0 0; 1 0 0; .5 .5 .5; [0.851 0.325 0.098]]; 
color = color(guseq,:);
close
figure1 = figure('Visible','off'); hold on;
for gi = guseq
    gn = regexprep(DataG(gi).name,'_',' ');
    d = DataG(gi).speedbm;
    y = nanmean(d);
    sd = nanstd(d);
    n = sum(~isnan(d));
    se = sd./sqrt(n-1);
    se2 = se.*2;
    x = DataG(gi).time(1,:);
    c = color(gi,:);
    e1 = plot(x,y,'Color',c,'Marker','none','DisplayName',gn,...
        'Linewidth',3);
end
xlim([min(x) max(x)])
legend1 = legend('show');
set(legend1,'Location','northeast','EdgeColor',[1 1 1]);
savename = sprintf('ephys t28-30');
printfig(savename,pSave,'w',4,'h',4,'closefig',1);


% find baseline
for gi = 1:size(DataG,2)
    % generate baseline
    i = (1:find(DataG(gi).time(1,:)==0)-1);
    d = DataG(gi).speedbm(:,i);
    DataG(gi).bs_data = d;
    [m,m_up,m_low,sd,se] = ephys_findbaseline(d);
    DataG(gi).bs_mean = m;
    DataG(gi).bs_sd = sd;
    DataG(gi).bs_se = se;
    DataG(gi).bs_upper = m_up;
    DataG(gi).bs_lower = m_low;
    fprintf('- baseline speed (%s): [%.3f]-[%.3f]\n',...
        char(DataG(gi).name),DataG(gi).bs_lower,DataG(gi).bs_upper);
end
% anova for baseline differences
a = nanmean(DataG(1).bs_data')';
b = nanmean(DataG(2).bs_data')';
x = [a;b];
group = [repmat({DataG(1).name},numel(a),1); repmat({DataG(2).name},numel(b),1)];
i = isnan(x);
x(i) = []; group(i) = [];
[text,T] = anova1_autoresults(x,group);
% store in stats
Stats(end+1).name = 'baseline';
Stats(end).anova = text;
Stats(end).descriptive = T;


% find mean rise peaks not on tap time(existance, value, time)
statsname = 'risepeak';
for gi = 1:size(DataG,2)
    % define range of time after tap
    bs = DataG(gi).bs_upper;
    time = DataG(gi).time(1,:);
    i = (find(time==0)+1):DataG(gi).timepts;
    d = nanmean(DataG(gi).speedbm(:,i)); % find mean
    t = time(:,i);
    [ymax,ymaxtime] = ephys_findrisepeak(d,t,bs);
    DataG(gi).(statsname) = ymax;
    DataG(gi).([statsname,'_time']) = ymaxtime;
end

