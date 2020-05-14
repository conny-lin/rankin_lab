function [str,v] = taptime_result(GraphData,ST,gname,statstext,varargin)

%% default setting
statxt = 'increased';
alpha = 0.05;
pvlimit = 0.001;
responseT1 = 0.1;
responseT2 = 1;
baselinetime = -0.1;

% VaraginProcesser


%%
xtime = [0.1:0.1:0.5];
rtime = xtime;

%%
fn = [gname,'_y'];
[rq,rd] = response_data2(GraphData,gname);
[pv,~] = get_pvalue(ST,10,gname);

switch statstext
    case 'increased'
        rval = rq.higher; 
        rtext = 'increased above baseline';
    case 'decreased'
        rval = rq.lower;
        rtext = 'decreased below baseline';

    case 'reversed'
        rval = rq.rev;
        rtext = 'decreased below zero';

end
        
[pstr,~,tv] = pvaluestring(rtime,rval,pv,alpha,pvlimit);
s = taptime_pvalue(tv,pstr);
if ~isempty(s)
    str = sprintf('significantly %s %s',rtext,s);
    v=true;
else
    statstext_mod = statstext(1:end-1);
    str = sprintf('did not significantly %s',rtext);
    v=false;

end













