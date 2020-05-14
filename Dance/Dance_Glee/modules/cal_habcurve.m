function [B,p] = cal_habcurve(D,pSave)
%% cal_habcurve
% input
%     D = MWTSet.Data.ByGroupPerPlate
%     pSave = save home folder

%% create save folder
pSaveA = [pSave,'/Stats HabCurve'];
if isdir(pSaveA) == 0; mkdir(pSaveA); end  


%% get data
gnames = fieldnames(D);
B = struct;
B.GroupNames = gnames;
% get measure list
msrlist = fieldnames(D.(gnames{1}));
% exclude ID, time and N measures
msrlist(...
            ismember(msrlist,{'MWTplateID','MWTind','time'}) | ...
            regexpcellout(msrlist,'SE\>') |...
            regexpcellout(msrlist,'^N')...
        ) = [];
% create output name
msrlistname = regexprep(msrlist,'_Mean','');
% create base entries (group name and plate N)
for g = 1:numel(gnames)
    gname = gnames{g};
    D1 = D.(gname);
    plateN = numel(D1.MWTplateID);
    B.PlateN(g,1) = plateN;
end


%% cal
for mi = 1:numel(msrlist)
    A = cell(1,size(gnames,1));
    T = table;
    for g = 1:numel(gnames)
        gname = gnames{g};
        D1 = D.(gname);
        X = (1:size(D1.N_TotalN,1))';
        msr = msrlist{mi};
        msrname = msrlistname{mi};
        d = D1.(msr);
        A{g} = d';
        B.(msrname).GroupName{g} = gname;
        B.(msrname).N{g} = size(d,2);
        B.(msrname).X(:,g) = X;
        B.(msrname).Y(:,g) = nanmean(d,2);
        if size(d,2) > 1
            sdn = nanstd(d')';
            sen = sdn./...
                sqrt(repmat(size(d,2),size(d,1),1));
        else
            sdn = nan(numel(X),1);
            sen = sdn;
        end
        B.(msrname).SD(:,g) = sdn;
        B.(msrname).E(:,g) = sen;
        % collect raw data
        t = table;
        t.mwtname = D.(gname).MWTplateID;
        t.groupname = cellstr(repmat(gname,size(t,1),1));
        a = regexp(gname,'_','split');
        a = a(:,1);
        t.strain = cellstr(repmat(a,size(t,1),1));
        a = regexp(gname,'(?<=^\w{1,}\d{1,}[_])\w*','match');
        if isempty(a) == 0
            t.dose = repmat(a,size(t,1),1);
        else
            t.dose = repmat({''},size(t,1),1);
        end
        t.data = cell(size(d,2),1);
        for x = 1:size(d,2)
            t.data{x} = d(:,x);
        end
        T = [T;t];
    end

    B.(msrname).Raw = T;
    % repeated measures anova
    if numel(A) > 1
        [p,t] = anova_rm(A,'off');
        B.(msrname).RMANOVA = t;
        b = anovan_textresult(t);
        B.(msrname).RMANOVA_textoutput = b;
        str = sprintf('%s/%s HabCurve RMANOVA.txt',pSaveA,msrname);
        fid = fopen(str,'w');
        fprintf(fid,'Repeated Measures ANOVA of %s all groups:\n',msrname);
        for x = 1:numel(b)
            fprintf(fid,'%s\n',b{x});
        end
        fclose(fid);
    end
    
    
    
    % run pair-wise rmanova if more than two groups
    if numel(gnames) > 2
        % construct comparison matrix
        n = numel(gnames);
        a = [];
        for x = 1:n-1
            v1 = repmat(x,n-x,1);
            v2 = (x+1:n)';
            b = [v1 v2];
            a = [a;b];
        end
        cmpindx = a;
        % rmanova
        b = {};
        for x = 1:size(cmpindx,1)
            a = cell(1,2);
            gn1 = gnames(cmpindx(x,1));
            gn2 = gnames(cmpindx(x,2));
            a{1} = cell2mat(T.data(ismember(T.groupname,gn1))');
            a{2} = cell2mat(T.data(ismember(T.groupname,gn2))');
            [p,t] = anova_rm(a,'off');
            if x == size(cmpindx,1)
                b = [b;...
                    [{sprintf('Hab curve %s vs %s',char(gn1),char(gn2))};
                    anovan_textresult(t)]];
            else
                b = [b;...
                    [{sprintf('Hab curve %s vs %s',char(gn1),char(gn2))};
                    anovan_textresult(t);{' '}]];
            end
        end
        B.(msrname).RMANOVA_PAIRWISE_textoutput = b;
        % save if pSave is provided
        if nargin == 2
           str = sprintf('%s/%s HabCurve RMANOVA pairwise.txt',pSaveA,msrname);
           fid = fopen(str,'w');
           fprintf(fid,'Repeated Measures ANOVA between groups:\n');
           for x = 1:numel(b)
               fprintf(fid,'%s\n',b{x});
           end
           fclose(fid);
        end
    end


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

