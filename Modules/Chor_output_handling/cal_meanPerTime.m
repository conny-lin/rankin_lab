function S = cal_meanPerTime(D,time,group,varargin)
%% instruction
% time has to be interval such as [0 10 20 30 40]. the output will
% has 4 time points: 10s = 0-10, 20s = 10-20, 30s = 20-30, 40s = 30-40



%% defaults
msr = D.Properties.VariableNames;
groupname = 'gname';
vararginProcessor;


%%
S = struct;
groupu = unique(group);
for mi = 1:numel(msr)
    TS = table;
    for gi = 1:numel(groupu)
        i = ismember(group,groupu(gi));
        T = grpstatsTable(D.(msr{mi})(i),time(i));
        i = ismember(T.Properties.VariableNames,'gnameu');
        T.Properties.VariableNames(i) = {'time'};
        T.time = cellfun(@str2num,T.time);     
        g = repmat(groupu(gi),size(T,1),1);
        T1 = table;
        T1.(groupname) = g;
        T1 = [T1 T];
        TS = [TS;T1];
    end
    S.(msr{mi}) = TS;
end


end













