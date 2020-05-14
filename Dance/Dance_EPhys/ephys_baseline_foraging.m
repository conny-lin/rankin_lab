% function ephys_baseline_foraging(BL)


%% examine sample
% examine sample size
[rn,cn] = size(BL);
fprintf('baseline foraging speed will be determined by:\n');
fprintf('- %d subjects x %d time points = %d data points\n',rn,cn,rn*cn);
if rn < 200
    fprintf('- subjects number < 200, need random sampling\n');
    randSample = true;
else
    randSample = false;
end


% examine reversals
revpct = (sum(sum(BL<0))./numel(BL)).*100;
fprintf('- reversal frequency: %.1f%% ',revpct);
revpct_acceptable = 10;
fprintf('(<%d%% ',revpct_acceptable);
if revpct < revpct_acceptable
   fprintf('acceptable)\n');
   exclude_reversal = true;
else
    fprintf('unacceptable)\n');
    exclude_reversal = false; 
    error('code to deal with this');
end


% examine pauses
pausepct = (sum(sum(BL==0))./numel(BL)).*100;
fprintf('- pause frequency: %.1f%% ',pausepct);
pausepct_acceptable = 20;
fprintf('(<%d%% ',pausepct_acceptable);

if pausepct < pausepct_acceptable
   fprintf('acceptable)\n');
   exclude_pause = true;
else
    fprintf('unacceptable)\n');
    exclude_pause = false; 
    error('code to deal with this');
end



%% exclude data
BLS = reshape(BL,numel(BL),1);
if exclude_reversal; BLS(BLS<0) = []; end
if exclude_pause; BLS(BLS==0) = []; end
BLS(isnan(BLS)) = []; 
n = numel(BLS);
fprintf('incuded data points: %d\n',n);


%% 2SD

sd = std(BLS);
se = sd/sqrt(n-1);
m = mean(BLS);
m1 = m-(sd);
m2 = m+(sd);
q1 = quantile(BLS,.25);
q2 = quantile(BLS,.75);
[q1 q2; m1 m2]
median(BLS)

















