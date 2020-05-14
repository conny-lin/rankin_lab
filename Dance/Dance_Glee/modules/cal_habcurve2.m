function [G,varcombo,B] = cal_habcurve2(D,varargin)
%% cal_habcurve2, create graph output and stats
% input
%     D = MWTSet.Data.trv2;
%     pSave = save home folder

%% DEFAULTS AND VARARGIN
var = {'groupname'}; % default group name
% process varargin
vararginProcessor; 


%% create var combintation
% check if input has var variables
vn = D.Properties.VariableNames;
if sum(ismember(var,vn)) ~= numel(var)
    error('some var input are not found');
else
    % make var table
    t = table;
    for x = 1:numel(var)
       t.(var{x}) = D.(var{x});
    end
    varref = t;
    varcombo = unique(t,'rows');
end


%% calculate habituation curve per time per var combo
msrlist = {'RevFreq','RevSpeed','RevDur'};
B = struct;
for mi = 1:numel(msrlist)
    msr =  msrlist{mi};
    A = table;
    for vi = 1:size(varcombo,1)
        % find group combo
        vc = varcombo(vi,:);
        i = ismember(varref,vc,'rows');
        % get data
        g = D.tap(i);
        X = D.(msr)(i);
        [gn,n,m,se] = grpstats(X,g,{'gname','numel','mean','sem'});
        d = [cellfun(@str2num,gn) n m se];
        d = array2table(d,'VariableNames',{'time','N','mean','SE'});
        % add legend
        leg = repmat(vc,size(d,1),1);
        d = [leg, d];
        A = [A;d];
    end
    B.(msr) = A;
end


%% make graphing variables (time as x)
% creae matrix where X = z dimension 1, Y = 2, E = 3, N = 4;
msrlist = {'RevFreq','RevSpeed','RevDur'};
empt = nan(numel(unique(D.tap)),size(varcombo,1),4);
G = struct;
for mi = 1:numel(msrlist)
    msr =  msrlist{mi};
    A = table;
    M = empt;
    for vi = 1:size(varcombo,1)
        % find group combo
        vc = varcombo(vi,:);
        i = ismember(varref,vc,'rows');
        % get data
        g = D.tap(i);
        d = D.(msr)(i);
        [gn,n,m,se] = grpstats(d,g,{'gname','numel','mean','sem'});
        M(:,vi,4) = n;
        M(:,vi,1) = cellfun(@str2num,gn);
        M(:,vi,2) = m;
        M(:,vi,3) = se;
    end
G.(msr) = M;
end






 
%% old code 
% 
%         return
%         %%
%         X = D1.RevFreq_Mean'
%         [p, table] = anova_rm(X,'on')
% 
% 
%         B.RevDur.X(:,g) = X;
%         B.RevDur.Y(:,g) = nanmean(D1.RevDur_Mean,2);
%         B.RevDur.SD(:,g) = nanstd(D1.RevDur_Mean')';
%         B.RevDur.E(:,g) = B.RevDur.SD(:,g)./...
%             sqrt(repmat(plateN,size(B.RevDur.SD(:,g),1),1));
% 
%         B.RevSpeed.X(:,g) = X;
%         B.RevSpeed.Y(:,g) = nanmean(D1.RevSpeed_Mean,2);
%         B.RevSpeed.SD(:,g) = nanstd(D1.RevSpeed_Mean')';
%         B.RevSpeed.E(:,g) = B.RevSpeed.SD(:,g)./...
%             sqrt(repmat(plateN,size(B.RevSpeed.SD(:,g),1),1));
%     end
% end

% for g = 1:numel(gnames)
%     gname = gnames{g};
% 
% 
%     D1 = D.(gname);
%     
%     
%     plateN = numel(D1.MWTplateID);
%     X = (1:size(D1.N_TotalN,1))';
% 
%     B.PlateN(g,1) = plateN;
%     for mi = 1:numel(msrlist)
%         msr = msrlist{mi};
%         msrname = msrlistname{mi};
%         
%         B.(msrname).X(:,g) = X;
%         B.(msrname).Y(:,g) = nanmean(D1.(msr),2);
%         B.(msrname).E(:,g) = nanstd(D1.(msr)')';
% 
%         B.N_Reversed.X(:,g) = X;
%         B.N_Reversed.Y(:,g) = nanmean(D1.N_Reversed,2);
%         B.N_Reversed.E(:,g) = nanmean(D1.N_Reversed,2);
% 
%         B.N_NoResponse.X(:,g) = X;
%         B.N_NoResponse.Y(:,g) = nanmean(D1.N_NoResponse,2);
%         B.N_NoResponse.E(:,g) = nanstd(D1.N_NoResponse')';
% 
%         B.RevFreq.X(:,g) = X;
%         B.RevFreq.Y(:,g) = nanmean(D1.RevFreq_Mean,2);
%         B.RevFreq.SD(:,g) = nanstd(D1.RevFreq_Mean')';
%         B.RevFreq.E(:,g) = B.RevFreq.SD(:,g)./...
%             sqrt(repmat(plateN,size(B.RevFreq.SD(:,g),1),1));
%         return
%         %%
%         X = D1.RevFreq_Mean'
%         [p, table] = anova_rm(X,'on')
% 
% 
%         B.RevDur.X(:,g) = X;
%         B.RevDur.Y(:,g) = nanmean(D1.RevDur_Mean,2);
%         B.RevDur.SD(:,g) = nanstd(D1.RevDur_Mean')';
%         B.RevDur.E(:,g) = B.RevDur.SD(:,g)./...
%             sqrt(repmat(plateN,size(B.RevDur.SD(:,g),1),1));
% 
%         B.RevSpeed.X(:,g) = X;
%         B.RevSpeed.Y(:,g) = nanmean(D1.RevSpeed_Mean,2);
%         B.RevSpeed.SD(:,g) = nanstd(D1.RevSpeed_Mean')';
%         B.RevSpeed.E(:,g) = B.RevSpeed.SD(:,g)./...
%             sqrt(repmat(plateN,size(B.RevSpeed.SD(:,g),1),1));
%     end
% end

