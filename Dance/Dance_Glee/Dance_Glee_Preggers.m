function MWTSet = Dance_Glee_Preggers(pMWT,varargin)
%% Dance_Glee_Preggers
% visualize individual plate by exp variations
% Dance_Glee: specifically designed to analyze alcohol effect on STH for N2 and
% N2_400mM
% 
% INPUTS:
%     pMWT = cell array of paths to MWT plates


%% DEFAULTS AND VARARGIN
pSave = '/Users/connylin/Dropbox/RL/Dance Output';
analysisNLimit = 20;
outputsuffix = '';
saveimport = 1;
%% STANDARD DANCE PROCESSING
% add general function paths
funPath = mfilename('fullpath');
funName = mfilename;
addpath([fileparts(fileparts(funPath)),'/Modules']); 
DanceM_funpath(funPath);

% process varargin
vararginProcessor; VarIn = varargin;
% create MWTSet
DanceM_MWTSetStd


%% [suspend] validate input: only one strain
% [~,gn] = cellfun(@fileparts,cellfun(@fileparts,pMWTS,'UniformOutput',0),'UniformOutput',0);
% a = unique(gn);
% a = regexpcellout(a,'_','split');
% strainU = unique(a(:,1));
% if numel(strainU) >1
%     error('This function can only accomodate one strain input');
% end
% refStrain = char(strainU);
% MWTSet.refStrain = refStrain;
% str = sprintf('(?<=^%s[_])\\w*',refStrain);
% condition = regexpcellout(gn,str,'match');
% conditionU = unique(condition(~cellfun(@isempty,condition)));
% if numel(conditionU) > 1; 
%     error('This function can only accomodate one condition input');
% end
% MWTSet.condition = conditionU;

%% CHOR
leg = chormaster4('ShaneSpark',pMWT);
MWTSet.chorlegend = leg;


%% TRV
% IMPORT .TRV (revised 20151126)
Data = import_trv(pMWT);
i = cellfun(@isempty,Data.data);
if sum(i) > 0
    Data_notrv = Data(i,:);
    export_MWTinfoTable(Data_notrv.mwtpath,'pSave',pSave,'suffix','no trv');
    Data(i,:) = [];
end

%% validation
A = struct;
[A.trv, A.trv_badtap] = validate_trv_tap(Data);
if isempty(A.trv_badtap) == 0
    export_MWTinfoTable(A.trv_badtap.mwtpath,'pSave',pSave,'suffix','bad tap');
end
[A.trv, A.trv_lowN] = validate_trv_N(A.trv,analysisNLimit);
if isempty(A.trv_lowN) == 0
    export_MWTinfoTable(A.trv_lowN.mwtpath,'pSave',pSave,'suffix',sprintf('N <%d',analysisNLimit));
end
[A.trv, A.trv_freqNaN] = validate_trv_freqNaN(A.trv);
if isempty(A.trv_freqNaN) == 0
    export_MWTinfoTable(A.trv_freqNaN.mwtpath,'pSave',pSave,'suffix','freq Nan');
end
export_MWTinfoTable(A.trv.mwtpath,'pSave',pSave,'suffix','included');
expsummaryTable(A.trv.mwtpath,'savetable',1,'pSave',pSave,'suffix','included','display',false);
if saveimport == 1; MWTSet.Import = A; end

%% transform trv
MWTSet.Data.trv = transform_trv_v0(A.trv);
% save raw data
export_trv_raw(pSave,MWTSet.Data.trv);


%% organize trv output 
% MWTSet.Data.trv_bygroup = orgdata_bygroup(MWTSet.Data.trv);

% by exp
% MWTSet.Data.trv_bygroup_byexp = orgdata_bygroup_byexp(MWTSet.Data.trv);


%% GRAPH: PER PLATE
pSaveA = [pSave,'/Graph Hab Curve by plates'];
if isdir(pSaveA) == 0; mkdir(pSaveA); end
D = MWTSet.Data.trv;
msru = fieldnames(D);
for mi = 1:numel(msru)
    msr = msru{mi};
    d = D.(msr);
    fnu = fieldnames(d);
    p = d.pMWT;
    DB = parseMWTinfo(p);
    gnu = unique(DB.groupname);
    figure1 = figure('Color',[1 1 1],'Visible','off');
    for gi = 1:numel(gnu)
        gn = gnu{gi};
        enu = unique(DB.expname(ismember(DB.groupname,gn)));
        for ei = 1:numel(enu)
            en = enu{ei};
            enName1 = char(regexp(en,'\<\d{8}[A-Z]{1}[_][A-Z]{2}[_]\d+s\d+x\d+s\d+s','match'));
            enName2 = ['E',enName1];
            i = ismember(DB.groupname,gn) & ismember(DB.expname,en);
            if sum(i) > 1
                % create col name
                % calculate hab curve
                d = D.(msr).mean(:,i);
                x = repmat((1:size(d,1))',1,size(d,2));

                if gi == 1
                    c = [0 0 0];
                    step = 0.5/(numel(enu)-1);
                    cstep = 0:step:0.5;
                    c = c + cstep(ei);
                elseif gi == 3
                    c = [0 0 1];
                    step = 0.8/(numel(enu)-1);
                    cstep = 0:step:0.8;
                    c(1) = c(1) + cstep(ei);
                    c(2) = c(2) + cstep(ei);
                elseif gi == 2
                    c = [1 0 0];
                    step = 1/(numel(enu)-1);
                    cstep = 0:step:1;
                    c(3) = c(3) + cstep(ei);
                end
                plot(x,d,'Color',c);
                hold on
            end
        end
    end
    ylabel(msr);
    xlabel('tap');
    savefigepsOnly150(msr,pSaveA);
end

%% CALCULATION: HAB CURVE (by group by exp)
% D = MWTSet.Data.trv;
% B = struct;
% msru = fieldnames(D);
% A = [];
% nRow = 1;
% G = table;
% for mi = 1:numel(msru)
%     msr = msru{mi};
%     d = D.(msr);
%     fnu = fieldnames(d);
%     p = d.pMWT;
%     DB = parseMWTinfo(p);
%     gnu = unique(DB.groupname);
%     enu = unique(DB.expname);
%     for gi = 1:numel(gnu)
%         for ei = 1:numel(enu)
%             gn = gnu{gi};
%             en = enu{ei};
%             enName1 = char(regexp(en,'\<\d{8}[A-Z]{1}[_][A-Z]{2}[_]\d+s\d+x\d+s\d+s','match'));
%             enName2 = ['E',enName1];
%             i = ismember(DB.groupname,gn) & ismember(DB.expname,en);
%             if sum(i) > 1
%                 % create col name
%                 % calculate hab curve
%                  d = D.(msr).mean(:,i)';
%                 A(nRow,:,1) = 1:size(d,2);
%                 A(nRow,:,2) = mean(d);
%                 A(nRow,:,3) = nanstd(d)./sqrt(sum(i)-1);
%                 
%                 % add group name
%                 group = table;
%                 group.groupname = {gn};
%                 group.expanme = {en};
%                 G = [G;group];
%                 
%                 nRow = nRow+1;
%             end
%         end
%     end
%     B.(msr).groupname = G;
%     B.(msr).graphdata = A;
% end
% MWTSet.Graph.HabCurve = B;
% 
% 
% % CALCULATION: HAB CURVE 
% % MWTSet.Graph.HabCurve = cal_habcurve(MWTSet.Data.ByGroupPerPlate,pSave);





%% SAVE MAT FILES
cd(pSave);
save([mfilename,'.mat'],'MWTSet');


%% Report done
fprintf('\n\n***DONE***\n\n');


end


%% [CODING] CALCULATION: ASYMPTOTIC LEVEL AS REACHING NO SIGN DIFF [CODING]
% this needs trinity data raw calculation
% 1 - asymptotic level is defined as no significant differences with the 
% % responses to next 3 taps)
% M = {'RevFreq', 'RevDur','RevSpeed'};
% pSaveA = [pSave,'/','Stats Hablevel stable'];
% if isdir(pSaveA) == 0; mkdir(pSaveA); end
% D = MWTSet.Data.ByGroupPerPlate;
% gnames = fieldnames(D);
% A = struct;
% 
% for m = 1:numel(M);% for each measure
%     % calculate hab level
%     B = []; G = {};
%     for g = 1:numel(gnames)
%         
%         data = D.(gnames{g}).([M{m},'_Mean']);
%         % get lowest response
% %         datadiff = data - repmat(min(data),size(data,1),1);
% %         % find area under the curve
% %         platearea = nan(size(datadiff,2),1);
%         for mwti = 1:size(data,2) 
%             % get trv data
%             data = MWTSet.Data.Import.trv{D.(gnames{g}).MWTind(mwti)};
%             d = data.([M{m}]);
%             e = data.([M{m},'_SE']);
%             
%             
%             return
%             a = datadiff(:,mwti);
%             dint = nan(size(a,2)-1,1);
%             for ti = 1:size(a,1)-1
%                 n1 = a(ti);
%                 n2 = a(ti+1);
%                 if n1 > n2; 
%                     nmax = n1; 
%                     nmin = n2;
%                 elseif n1==n2; 
%                     nmax = n1;
%                     nmin = n2;
%                 else
%                     nmax = n2;
%                     nmin= n1;
%                 end
%                 basearea = nmin*1;
%                 trianglearea = (nmax*1)/2;
%                 totalarea = basearea + trianglearea;
%                 dint(ti) = totalarea;      
%             end
%             platearea(mwti) = sum(dint);
%         end
%                 
%         % get data
%         B = [B;platearea];
%         n = numel(platearea);
%         G = [G; repmat(gnames(g), n,1)];
%     end
%     
%     
%     % stats
%     [m2, n2, se2,gnames2] = grpstats(B,G,{'mean','numel','sem','gname'});
%     A.(M{m}).GroupName = gnames2;
%     A.(M{m}).N = n2;
%     A.(M{m}).Y = m2;
%     A.(M{m}).E = se2; 
%     if numel(unique(G)) > 1
%         % anova
%         [p,t,stats] = anova1(B,G,'off');
%         [c,m1,h,gnames] = multcompare(stats,'ctype','bonferroni','display','off');
%         A.(M{m}).ANOVA = t;
%         % export anova
%         vname = regexprep(t(1,:),'Prob>F','P_value');
%         T = cell2table(t(2:end,:),'VariableNames',vname);
%         cd(pSaveA); 
%         writetable(T,[M{m},' Initial ANOVA.csv'],'Delimiter',',');
%         
%         % posthoc
%         i = ((c(:,2) >0 & c(:,4) >0) + (c(:,2) <0 & c(:,4) <0)) >0;
%         a = [gnames(c(:,1)), gnames(c(:,2))];
%         a(i,3) = {'< 0.05'};
%         a(~i,3) = {'not significant'};
%         A.(M{m}).posthoc = a;
%         % export posthoc
%         T = cell2table(a);
%         cd(pSaveA); 
%         writetable(T,[M{m},' Initial posthoc bonferroni.csv'],'Delimiter',',');
%         
%     end
% end
% MWTSet.Graph.HabRate_integral = A;

% Habituation rate (how fast does worms habituated to stable asymptote level; 

%% GRAPH: HABITUATION CURVES
% pSaveA = [pSave,'/Graph Habituation curves'];
% if isdir(pSaveA) == 0; mkdir(pSaveA); end
% M = {'N_Total'
%     'N_Reversed'
%     'N_NoResponse'
%     'RevFreq'
%     'RevDur'
%     'RevSpeed'};
% Graph = MWTSet.Graph.HabCurve;
% GroupName = regexprep(Graph.GroupNames,'_',' ');
% timestamp = MWTSet.timestamp;
% % get N2 
% N2i = ~cellfun(@isempty,regexp(GroupName,'^N2'));
% 
% for m = 1:numel(M);% for each measure
%     % get coordinates
%     X = Graph.(M{m}).X;
%     Y = Graph.(M{m}).Y;
%     E = Graph.(M{m}).E;
% 
%     % create output table legend
%     a = cellstr([repmat('t',size(Y,1),1) num2str([1:size(Y,1)]')]);
%     a = regexprep(a,' ','');
%     b = cellstr([repmat('SE',size(Y,1),1) num2str([1:size(Y,1)]')]);
%     b = regexprep(b,' ','');
%     vnames = [{'groupname'}; a;b];   
%     % make table 
%     d = [GroupName num2cell(Y') num2cell(E')];
%     T = cell2table(d,'VariableNames',vnames);
%     cd(pSaveA);
%     writetable(T,[M{m},'.csv'],'Delimiter',',');
%         
%     % Create figure
%     figure1 = figure('Color',[1 1 1],'Visible','off');
%     axes1 = axes('Parent',figure1,'FontSize',18);
%     box(axes1,'off');
%     xlim(axes1,[0 size(Y,1)+1]);
%     ylim(axes1,[min(min(Y-E))*.9 max(max(Y+E))*1.1]);
%     hold(axes1,'all');
%     errorbar1 = errorbar(X,Y,E,'LineWidth',1);
%     color = [0 0 0; 1 0 0; [0.5 0.5 0.5]; [0.04 0.52 0.78]]; 
%     n = min([size(color,1) size(Y,2)]);
%     for x = 1:n
%         set(errorbar1(x),'Color',color(x,1:3))
%     end
% 
%     for x = 1:size(Y,2)
%         set(errorbar1(x),'DisplayName',GroupName{x})
%     end
%     titlestr = timestamp;
%     title(titlestr,'FontSize',12);
%     xlabel('Tap','FontSize',18);
%     ylabel(regexprep(M{m},'_',' '),'FontSize',18);
%     legend1 = legend(axes1,'show');
%     set(legend1,'EdgeColor',[1 1 1],...
%         'Location','NorthEastOutside',...
%         'YColor',[1 1 1],'XColor',[1 1 1],'FontSize',12);
% 
%     % text strings
%     a = Graph.PlateN;
%     b = '';
%     for x = 1:numel(a)
%         b = [b,num2str(a(x)),', '];
%     end
%     str1 = ['N(plate) = ',b(:,1:end-2)];
%     textboxstr = [str1];
%     annotation(figure1,'textbox',...
%         [0.70 0.015 0.256 0.05],...
%         'String',textboxstr,...
%         'FitBoxToText','off',...
%         'EdgeColor',[1 1 1]);
%     % save fig
%     savefigepsOnly(M{m},pSaveA);
%     
% end












