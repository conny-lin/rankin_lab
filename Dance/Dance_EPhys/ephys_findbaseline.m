function [m,m_up,m_low,sd,se] = ephys_findbaseline(d)

% for gi = 1:size(DataG,2)
    % generate baseline
%     i = (1:find(DataG(gi).time(1,:)==0)-1);
%     DataG(gi).bs_data = DataG(gi).speedbm(:,i);
m = nanmean(nanmean(d));
sd = nanstd(nanmean(d));
se = nanstd(nanmean(d))./sqrt(numel(nanmean(d))-1);
m_up = m+(sd*2);
m_low = m-(sd*2);
% fprintf('- baseline speed (%s): [%.3f]-[%.3f]\n',char(DataG(gi).name),DataG(gi).bs_lower,DataG(gi).bs_upper);
% end