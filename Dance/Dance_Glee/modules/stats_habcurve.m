function [S] = stats_habcurve(D,VInd,pSave, varargin)
%% stats_habcurve
% input
%     D = MWTSet.Data.ByGroupPerPlate
%     pSave = save home folder


%% DEFAULTS AND VARARGIN
var = {'groupname'}; % default group name
% var = {'strain','rx'}; 
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

%% make literal legend
leg = DanceM_convert_MWTInd2text(varcombo,VInd);

%% make graphing variables (time as x)
% creae matrix where X = z dimension 1, Y = 2, E = 3, N = 4;
msrlist = {'RevFreq','RevSpeed','RevDur'};
% empt = nan(numel(unique(D.tap)),size(varcombo,1),4);
S = struct;
str = sprintf('%s/HabCurve RMANOVA.txt',pSave);
fid = fopen(str,'w');
fprintf(fid,'Repeated Measures ANOVA for %s*time:\n',strjoin(var,'*'));

for mi = 1:numel(msrlist)
    msr =  msrlist{mi};
    A = {};
    val = true(size(varcombo,1),1);
    B = {};
    % make literal legend
    leg = DanceM_convert_MWTInd2text(varcombo,VInd);
    for vi = 1:size(varcombo,1)
        % find group combo
        vc = varcombo(vi,:);
        i = ismember(varref,vc,'rows');
        % get data
        t = D.tap(i);
        d = D.(msr)(i);
        % get plate
        pt = D.mwtname(i);
        ptu = unique(pt);
        a = nan(numel(ptu),numel(unique(t)));
        for pti = 1:numel(ptu)
            j = pt == ptu(pti);
            col = t(j);
            b = d(j);
            a(pti,col) = b;
        end
        B{vi} = a;

        % get rid of nan data
        a(any(isnan(a)'),:) = [];
        if size(a,1) <2
            val(vi) = false;
            a = [];
        end
        A{vi} = a;
    end
    A(~val) = [];
    
    if sum(val) ~= size(varcombo,1)
      u = table2array(leg(~val,:));
      newleg = cell(size(u,2),1);
      for x = 1:size(u,1)
         newleg{x} = strjoin(u(x,:),'*');
      end
    end

    % repeated measures anova
    if numel(A) > 1
        [p,t] = anova_rm(A,'off');
        S.(msr).RMANOVA = t;
        b = anovan_textresult(t);
        fprintf(fid,'%s ------\n',msr);
        for x = 1:numel(b)
            fprintf(fid,'%s\n',b{x});
        end
        
        if sum(val) ~= size(varcombo,1)
            fprintf(fid,'*excluded groups: ');
           for x = 1:numel(newleg)
                fprintf(fid,'%s ',newleg{x});
           end 
            fprintf(fid,'\n');
            
        end
        fprintf(fid,'\n');
    end
end
fclose(fid);



% %% pair wise
%     
%     
%     % run pair-wise rmanova if more than two groups
%     if numel(gnames) > 2
%         % construct comparison matrix
%         n = numel(gnames);
%         a = [];
%         for x = 1:n-1
%             v1 = repmat(x,n-x,1);
%             v2 = (x+1:n)';
%             b = [v1 v2];
%             a = [a;b];
%         end
%         cmpindx = a;
%         % rmanova
%         b = {};
%         for x = 1:size(cmpindx,1)
%             a = cell(1,2);
%             gn1 = gnames(cmpindx(x,1));
%             gn2 = gnames(cmpindx(x,2));
%             a{1} = cell2mat(T.data(ismember(T.groupname,gn1))');
%             a{2} = cell2mat(T.data(ismember(T.groupname,gn2))');
%             [p,t] = anova_rm(a,'off');
%             if x == size(cmpindx,1)
%                 b = [b;...
%                     [{sprintf('Hab curve %s vs %s',char(gn1),char(gn2))};
%                     anovan_textresult(t)]];
%             else
%                 b = [b;...
%                     [{sprintf('Hab curve %s vs %s',char(gn1),char(gn2))};
%                     anovan_textresult(t);{' '}]];
%             end
%         end
%         B.(msrname).RMANOVA_PAIRWISE_textoutput = b;
%         % save if pSave is provided
%         if nargin == 2
%            str = sprintf('%s/%s HabCurve RMANOVA pairwise.txt',pSaveA,msrname);
%            fid = fopen(str,'w');
%            fprintf(fid,'Repeated Measures ANOVA between groups:\n');
%            for x = 1:numel(b)
%                fprintf(fid,'%s\n',b{x});
%            end
%            fclose(fid);
%         end
%     end
% 
% 
% end

 
