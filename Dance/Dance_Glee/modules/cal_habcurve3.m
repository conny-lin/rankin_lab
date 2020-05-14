function M = cal_habcurve3(X,Y,G)
%% cal_habcurve2, create graph output and stats
% input
%     D = MWTSet.Data.trv2;
%     pSave = save home folder

% msrlist = {'RevFreq','RevSpeed','RevDur'};


%% calculate habituation curve per time per var combo (suspend)
% B = struct;
% for mi = 1:numel(msrlist)
%     msr =  msrlist{mi};
%     A = table;
%     for vi = 1:size(varcombo,1)
%         % find group combo
%         vc = varcombo(vi,:);
%         i = ismember(varref,vc,'rows');
%         
%         % get data
%         g = D.tap(i);
%         X = D.(msr)(i);
%         
%         [gn,n,m,se] = grpstats(X,g,{'gname','numel','mean','sem'});
%         d = [cellfun(@str2num,gn) n m se];
%         d = array2table(d,'VariableNames',{'time','N','mean','SE'});
%         
%         % add legend
%         leg = repmat(vc,size(d,1),1);
%         d = [leg, d];
%         A = [A;d];
%     end
%     B.(msr) = A;
% end


%% make graphing variables (time as x)
% creae matrix where X = z dimension 1, Y = 2, E = 3, N = 4;
% msrlist = {'RevFreq','RevSpeed','RevDur'};
GU = unique(G);
empt = nan(numel(unique(X)),numel(GU),4);
% G = struct;
% for mi = 1:numel(msrlist)
%     msr =  msrlist{mi};
%     A = table;
M = empt;
for gi = 1:numel(GU)
    % find group combo
    gname = GU(gi);
    i = ismember(G,gname);

    % get data
    g = X(i);
    d = Y(i);

    [gn,n,m,se] = grpstats(d,g,{'gname','numel','mean','sem'});
    M(:,gi,4) = n;
    M(:,gi,1) = cellfun(@str2num,gn);
    M(:,gi,2) = m;
    M(:,gi,3) = se;
end





 

