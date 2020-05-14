function MWTSet = Dance_Glee_VitaminD(pMWT,varargin)
%% INFORMATION
% specifically designed to analyze alcohol effect on STH for N2 and
% N2_400mM only
% comparison between a reference strain and a test strain
% statistics = per group by plates 
% organize output by groups
% pMWTS = path to MWT folders, output will be organized by group folders
% 
% optional inputs
%     pSave
%     analysisNLimit = 5; put minimum number of worms must be reversing to
%     first tap in order to qualify for analysis
%     alcoholTestDose
%     refStrain = 'N2';

%% DEFAULTS AND VARARGIN
pSave = '/Users/connylin/Dropbox/RL/Dance Output';
analysisNLimit = 5;
outputsuffix = '';
saveimport = 0;
GroupVar = {'groupname'}; % default group name
chorName = 'ShaneSpark';
MsrList = {'RevFreq','RevSpeed','RevDur'}; % define measure list 
%% STANDARD DANCE PROCESSING
% add general function paths
funPath = mfilename('fullpath');
funName = mfilename;
addpath([fileparts(fileparts(funPath)),'/Modules']); 
DanceM_funpath(funPath);
% process varargin
vararginProcessor; 
% create MWTSet
MWTSet = DanceM_MWTSetStd2(pMWT,funName,GroupVar,varargin,pSave,outputsuffix);

%
% end of section
%% CHOR
leg = chormaster4(chorName,pMWT);
% store legend in MWTSet
MWTSet.chorlegend = leg;

% end of section
%% IMPORT TRV
% get variables
MWTDBInd = MWTSet.Info.MWTDBInd;
VInd = MWTSet.Info.VInd;
GroupVarName = MWTSet.Info.GroupVarName;

% - import trv
Data = import_trv(pMWT);
i = cellfun(@isempty,Data.data);
if sum(i) > 0
    Data_notrv = Data(i,:);
    export_MWTinfoTable(Data_notrv.mwtpath,'pSave',pSave,'suffix','no trv');
    Data(i,:) = [];
end

% - trv validation
A = struct;

% validate tap N
[A.trv, A.trv_badtap] = validate_trv_tap(Data);
if isempty(A.trv_badtap) == 0
    export_MWTinfoTable(A.trv_badtap.mwtpath,'pSave',pSave,'suffix','bad tap');
end

% exclude plates with <analysisNlimit to first tap
[A.trv, A.trv_lowN] = validate_trv_N(A.trv,analysisNLimit);
if isempty(A.trv_lowN) == 0
    export_MWTinfoTable(A.trv_lowN.mwtpath,'pSave',pSave,'suffix',sprintf('N <%d',analysisNLimit));
end

% exclude plates with any missing frequency data
[A.trv, A.trv_freqNaN] = validate_trv_freqNaN(A.trv);
if isempty(A.trv_freqNaN) == 0
    export_MWTinfoTable(A.trv_freqNaN.mwtpath,'pSave',pSave,'suffix','freq Nan');
end

% export data information to csv
export_MWTinfoTable(A.trv.mwtpath,'pSave',pSave,'suffix','included');
expsummaryTable(A.trv.mwtpath,'savetable',1,'pSave',pSave,'suffix','included','display',false);
% store results in MWTSet
if saveimport == 1; MWTSet.Import = A; end
TRV = A.trv;

%% transform trv
A = A.trv;
VInd = MWTSet.Info.VInd;
[D,PF] = transform_trv_v2(A,'tapindx',1);
t = table;
[~,i] = ismember(PF,VInd.mwtpath);
t.mwtpath = i;
[~,j] = ismember(i,MWTDBInd.mwtpath);
t.(GroupVarName) = j;
D = [t D];
MWTSet.Data.trv2 = D;
% save raw data 
B = DanceM_convert_MWTInd2text(D,VInd); % conver index to text
writetable(B,sprintf('%s/Data trv.csv',pSave));

% end of section



%% CALCULATION: HAB CURVE


%% make graphing variables (time as x)
% creae matrix where X = z dimension 1, Y = 2, E = 3, N = 4;
% msrlist = {'RevFreq','RevSpeed','RevDur'};
D = MWTSet.Data.trv2;
A = struct;
for mi = 1:numel(MsrList)
    msr =  MsrList{mi};
    A = table;
    X = D.tap;
    Y = D.(msr);
    G = D.(GroupVarName);
    A.(msr) = cal_habcurve3(X,Y,G);
end
MWTSet.Graph.HabCurve = A;
return


%% START CODING HERE *****************************************************
% SAVE AND FINISH 
cd(pSave);
save([mfilename,'.mat'],'MWTSet');
% Report done
fprintf('\n\n***DONE***\n\n');
return


%% calculate stats: hab curve
D = MWTSet.Data.trv2;
VInd = MWTSet.Info.VInd;
%% calculate stats: hab curve
D = MWTSet.Data.trv2;
VInd = MWTSet.Info.VInd;
%% make graphing variables (time as x)
% creae matrix where X = z dimension 1, Y = 2, E = 3, N = 4;
MsrList = {'RevFreq','RevSpeed','RevDur'};
% empt = nan(numel(unique(D.tap)),size(varcombo,1),4);
S = struct;
str = sprintf('%s/HabCurve RMANOVA.txt',pSave);
fid = fopen(str,'w');
fprintf(fid,'Repeated Measures ANOVA for %s*time:\n',strjoin(VarNames,'*'));

for mi = 1:numel(MsrList)
    msr =  MsrList{mi};
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

%% CALCULATION: INITIAL
pSaveA = [pSave,'/','Stats Initial'];
prefix = 'Initial x strain x alcohol';
M = {'RevFreq', 'RevDur', 'RevSpeed'};
if isdir(pSaveA) == 0; mkdir(pSaveA); end
D = MWTSet.Data.ByGroupPerPlate;
gnames = fieldnames(D);
A = struct;
for m = 1:numel(M);% for each measure
    B = []; G = {}; mwtnames = {};
    for g = 1:numel(gnames)
        d1 = D.(gnames{g}).([M{m},'_Mean'])(1,:);
        B = [B;d1'];
        n = numel(d1);
        G = [G; repmat(gnames(g), n,1)];
        mwtnames = [mwtnames;D.(gnames{g}).MWTplateID];
    end
    % exploratory stats-------------------------------------------------
    [m2, n2, se2,gnames2] = grpstats(B,G,{'mean','numel','sem','gname'});
    A.(M{m}).GroupName = gnames2;
    A.(M{m}).N = n2;
    A.(M{m}).Y = m2;
    A.(M{m}).E = se2; 
    % stat test---------------------------------------------------------
    A = statsoutput_anova_alcohol(prefix,A,B,G,mwtnames,pSaveA,M{m});    
end
MWTSet.Graph.Initial = A;
%% CALCULATION: HAB LEVEL LAST 3 TAPS 
pSaveA = [pSave,'/','Stats HabLevel'];
prefix = 'Hablevel x strain x alcohol';
M = {'RevFreq', 'RevDur', 'RevSpeed'};
if isdir(pSaveA) == 0; mkdir(pSaveA); end
D = MWTSet.Data.ByGroupPerPlate;
gnames = fieldnames(D);
A = struct;
for m = 1:numel(M);% for each measure
    B = []; G = {}; mwtnames = {};
    for g = 1:numel(gnames)
        d = D.(gnames{g}).([M{m},'_Mean'])(end-2:end,:);
        % calculate last 3 taps average per plate
        d1 = nanmean(d);
        % get data
        B = [B;d1'];
        n = numel(d1);
        G = [G; repmat(gnames(g), n,1)];
        mwtnames = [mwtnames;D.(gnames{g}).MWTplateID];
    end
    % exploratory stats-------------------------------------------------
    [m2, n2, se2,gnames2] = grpstats(B,G,{'mean','numel','sem','gname'});
    A.(M{m}).GroupName = gnames2;
    A.(M{m}).N = n2;
    A.(M{m}).Y = m2;
    A.(M{m}).E = se2; 
    % stat test---------------------------------------------------------
    A = statsoutput_anova_alcohol(prefix,A,B,G,mwtnames,pSaveA,M{m});    
end
MWTSet.Graph.HabLevel = A;
%% CALCULATION HAB RATE: Andrew's rate - half life
% = first tap achieved half of hab level 
% (hab level (mean of last 3 taps) - initial)/2
pSaveA = [pSave,'/','Stats Habrate half life'];
prefix = 'HabRate_halflife x strain x alcohol';
M = {'RevFreq', 'RevDur', 'RevSpeed'};
if isdir(pSaveA) == 0; mkdir(pSaveA); end
D = MWTSet.Data.ByGroupPerPlate;
gnames = fieldnames(D);
A = struct;
for m = 1:numel(M);% for each measure
    B = []; G = {}; mwtnames = {};
    for g = 1:numel(gnames)
        data = D.(gnames{g}).([M{m},'_Mean']);
        % get initial
        init = data(1,:);
        % get hab level
        hab = nanmean(data(end-2:end,:)); 
        a = data - repmat(hab,size(data,1),1); % minus hab level
        b = a./repmat(a(1,:),size(a,1),1); % divided by initial
        i = b <= .5; % find 50% point
        midtap = nan(size(i,2),1);
        for x = 1:size(i,2)
            a = find(i(:,x));
            if isempty(a) == 0
                midtap(x) =  a(1);
            end
        end
        % get data
        B = [B;midtap];
        n = numel(midtap);
        G = [G; repmat(gnames(g), n,1)];
        mwtnames = [mwtnames;D.(gnames{g}).MWTplateID];
    end
    % exploratory stats-------------------------------------------------
    [m2, n2, se2,gnames2] = grpstats(B,G,{'mean','numel','sem','gname'});
    A.(M{m}).GroupName = gnames2;
    A.(M{m}).N = n2;
    A.(M{m}).Y = m2;
    A.(M{m}).E = se2; 
    % stat test---------------------------------------------------------
    A = statsoutput_anova_alcohol(prefix,A,B,G,mwtnames,pSaveA,M{m});    
end
MWTSet.Graph.HabRate_halflife = A;
%% CALCULATION HAB RATE: AREA BELOW THE CURVE (INTEGRAL)
% makes no assumption on which is the habituated level
% find the lowest point of response and set it as zero
% calculate area under the curve
% the smaller the area, the faster the animal reaches lowest level
% = first tap achieved half of hab level 
% (hab level (mean of last 3 taps) - initial)/2
pSaveA = [pSave,'/','Stats Habrate integral'];
prefix = 'HabRate_integral x strain x alcohol';
M = {'RevFreq', 'RevDur', 'RevSpeed'};
if isdir(pSaveA) == 0; mkdir(pSaveA); end
D = MWTSet.Data.ByGroupPerPlate;
gnames = fieldnames(D);
A = struct;
for m = 1:numel(M);% for each measure
    B = []; G = {}; mwtnames = {};
    for g = 1:numel(gnames)
        mwtnames = [mwtnames;D.(gnames{g}).MWTplateID];
        data = D.(gnames{g}).([M{m},'_Mean']);
        % get lowest response
        datadiff = data - repmat(min(data),size(data,1),1);
        % find area under the curve
        platearea = nan(size(datadiff,2),1);
        for mwti = 1:size(datadiff,2) 
            a = datadiff(:,mwti);
            dint = nan(size(a,2)-1,1);
            for ti = 1:size(a,1)-1
                n1 = a(ti);
                n2 = a(ti+1);
                if n1 > n2; 
                    nmax = n1; 
                    nmin = n2;
                elseif n1==n2; 
                    nmax = n1;
                    nmin = n2;
                else
                    nmax = n2;
                    nmin= n1;
                end
                basearea = nmin*1;
                trianglearea = (nmax*1)/2;
                totalarea = basearea + trianglearea;
                dint(ti) = totalarea;      
            end
            platearea(mwti) = sum(dint);
        end
        % get data
        B = [B;platearea];
        n = numel(platearea);
        G = [G; repmat(gnames(g), n,1)];
    end
    % exploratory stats-------------------------------------------------
    [m2, n2, se2,gnames2] = grpstats(B,G,{'mean','numel','sem','gname'});
    A.(M{m}).GroupName = gnames2;
    A.(M{m}).N = n2;
    A.(M{m}).Y = m2;
    A.(M{m}).E = se2; 
    % stat test---------------------------------------------------------
    A = statsoutput_anova_alcohol(prefix,A,B,G,mwtnames,pSaveA,M{m});    
end

MWTSet.Graph.HabRate_integral = A;
%% [suspend] STATS SUMMARY: strain vs strain_condition only (r20151201)
% % rownames
% pSaveA = pSave;
% % strainlist = [{refStrain};{testStrain}];
% M = {'RevDur','RevFreq','RevSpeed'};
% Assays = fieldnames(MWTSet.Graph);
% Assays(ismember(Assays,'HabCurve')) = [];
% T = table;
% a = {};
% for mi = 1:numel(M);
%     a = [a;repmat(M(mi),numel(Assays),1)];
% end
% T.msr = a;
% T.assays = repmat(Assays,numel(M),1);
% T.condition = cell(size(T,1),1);
% % T.(strainlist{2}) = cell(size(T,1),1);
% % T.ref_alcoholeffect_vs_ref = cell(size(T,1),1);
% % T.test_vs_ref_0mM = cell(size(T,1),1);
% 
% for ai = 1:numel(Assays)
%     for mi = 1:numel(M)
%         asr = Assays{ai};
%         msr = M{mi};
%         % repeat for test and ref strain
%         t = MWTSet.Graph.(asr).(msr).posthoc;
%         
%         return
%         
%         %%------------ need to fix
%         for si = 1:numel(strainlist)
%             i = ismember(t.group1,strainlist{si}) &...
%                 ismember(t.group2,sprintf('%s_%dmM',strainlist{si},alcoholTestDose));
%             if sum(i) == 0;
%                 T.(strainlist{si})(ismember(T.msr,msr) & ismember(T.assays,asr)) = {'n.s.'};
%             end
%             if sum(i) == 1
%                 n1 = MWTSet.Graph.(asr).(msr).Y(ismember(MWTSet.Graph.(asr).(msr).GroupName,strainlist{si}));
%                 n2 = MWTSet.Graph.(asr).(msr).Y(ismember(MWTSet.Graph.(asr).(msr).GroupName,sprintf('%s_%dmM',strainlist{si},alcoholTestDose)));
%                 a = n2-n1;
%                 if a < 0
%                     T.(strainlist{si})(ismember(T.msr,msr) & ismember(T.assays,asr)) = {'decrease'};
%                 elseif a> 0
%                     T.(strainlist{si})(ismember(T.msr,msr) & ismember(T.assays,asr)) = {'increase'};
%                 end
%             end
%         end
%         
%         %% ref_vs_test_0mM
%         name1 = refStrain; 
%         name2 = testStrain;
%         ti = ismember(t.group1,name1) & ismember(t.group2,name2) | ...
%             ismember(t.group1,name2) & ismember(t.group2,name1);
%         if sum(ti) == 0;
%             T.test_vs_ref_0mM(ismember(T.msr,msr) & ismember(T.assays,asr)) = {'n.s.'};
%         elseif sum(ti) == 1
%             n1 = MWTSet.Graph.(asr).(msr).Y(ismember(MWTSet.Graph.(asr).(msr).GroupName,name1));
%             n2 = MWTSet.Graph.(asr).(msr).Y(ismember(MWTSet.Graph.(asr).(msr).GroupName,name2));
%             a = n2-n1;
%             if a < 0
%                 T.test_vs_ref_0mM(ismember(T.msr,msr) & ismember(T.assays,asr)) = {'decrease'};
%             elseif a> 0
%                 T.test_vs_ref_0mM(ismember(T.msr,msr) & ismember(T.assays,asr)) = {'increase'};
%             end
%         end
%         
%     end
% end
% 
% % do alcohol effect analysis
% a = T.ref_alcoholeffect_vs_ref;
% for x = 1:size(T,1)
%    cond_ctrl = T.(refStrain)(x);
%    cond_test = T.(testStrain)(x);
%    if strcmp(cond_ctrl,'n.s.') == 1 
%        if strcmp(cond_test,'n.s.') == 1
%            a{x} = '--';
%        else
%            a{x} = 'sensitize';
%        end
%    else
%       if strcmp(cond_test,'n.s.') == 1
%           a{x} = 'ko';
%       elseif strcmp(cond_test,cond_ctrl) == 0
%           a{x} = 'opposite';
%       end
%    end
% end
% T.ref_alcoholeffect_vs_ref = a;
% MWTSet.AlcoholSummary = T;
% writetable(T,[pSaveA,'/AlcoholEffectSummary.csv']);




return
%% GRAPH: HABITUATION CURVES BY STRAINS + STATS
pSaveA = [pSave,'/Graph HabCurve'];
if isdir(pSaveA) == 0; mkdir(pSaveA); end
M = {'RevFreq', 'RevDur','RevSpeed'};
Graph = MWTSet.Graph.HabCurve;
HabLevel = MWTSet.Graph.HabLevel;
Initial = MWTSet.Graph.Initial;
HabCurve = MWTSet.Graph.HabCurve;
HabRate = MWTSet.Graph.HabRate_integral;
GroupName = regexprep(Graph.GroupNames,'_',' ');
groupN = numel(unique(GroupName));

N2i = ismember(GroupName,refStrain); % get N2 
color = [0 0 0; 1 0 0; [0.04 0.52 0.78]; [0.47843137383461 0.062745101749897 0.894117653369904]]; 
strain = refStrain;

for m = 1:numel(M);% for each measure
    msr = M{m};
    % get coordinates
    X = Graph.(M{m}).X;
    Y = Graph.(M{m}).Y;
    E = Graph.(M{m}).E;
    % make figures by group and N2
    X1 = X; Y1 = Y;  E1 = E;
  
    %% create habituation curve
    iG = ~N2i;
    idata = [find(N2i);find(iG)];
    gname = GroupName([find(N2i);find(iG)]);
    X = repmat(X1(:,1),1,sum(N2i)+sum(iG)); 
    Y = Y1(:,idata); 
    E = E1(:,idata);
    tapN = size(X,1);
    figure1 = figure('Color',[1 1 1],'Visible','off');
    axes1 = axes('Parent',figure1,'FontSize',18);
    box(axes1,'off');
    % x lim
    tapXlim = max(max(X))+1;
    xlim(axes1,[0 tapXlim+6]);
    % y lim
    if isnan(max(max(E))) == 1
        rspMin = min(min(Y))*.9;
        rspMax = max(max(Y));
    else
        rspMin = min(min(Y-E))*.9;
        rspMax = max(max(Y+E));
    end
    if strcmp(M{m},'RevSpeed') == 1
        ylim(axes1,[rspMin rspMax*1.1]);
    elseif strcmp(M{m},'RevSpeed') == 0
        ylim(axes1,[0 rspMax*1.1]);
    end
    hold(axes1,'all');
    errorbar1 = errorbar(X,Y,E,'LineWidth',1,'Marker','o');
    n = min([size(color,1) size(Y,2)]);
    for x = 1:n
        set(errorbar1(x),'Color',color(x,1:3),...
            'MarkerFaceColor',color(x,1:3),...
        'MarkerEdgeColor',color(x,1:3),...
        'LineStyle','none')
    end
    for x = 1:size(Y,2)
        set(errorbar1(x),'DisplayName',gname{x})
    end

    %% text strings
    a = Graph.PlateN(idata);
    b = '';
    for x = 1:numel(a)
        b = [b,num2str(a(x)),', '];
    end
    str1 = ['N(plate) = ',b(:,1:end-2)];
    title(str1,'FontSize',8);
    xlabel('Tap','FontSize',18);
    ylabel(regexprep(M{m},'_',' '),'FontSize',24);
    legend1 = legend(axes1,'show');
    set(legend1,'EdgeColor',[1 1 1],...
        'Location','NorthEastOutside',...
        'YColor',[1 1 1],'XColor',[1 1 1],'FontSize',12);


    %% graph initial dot
%             iG = ~cellfun(@isempty,regexp(MWTSet.Graph.Initial.(M{m}).GroupName,['^',testStrain]));
    N2i = ~cellfun(@isempty,regexp(MWTSet.Graph.Initial.(M{m}).GroupName,['^',refStrain])); % get N2
    % graph N2
    Y = []; E = []; X = [];
    Y(:,1) = Initial.(M{m}).Y(N2i);
    E(:,1) = Initial.(M{m}).E(N2i);
    X = repmat(tapXlim+1,size(Y));
    for x = 1:numel(Y)
        errorbar(X(x),Y(x),E(x),'MarkerSize',6,'Marker','o',...
        'LineStyle','none','LineWidth',1,'Color',color(x,1:3),...
        'MarkerFaceColor','none',...
        'MarkerEdgeColor',color(x,1:3));
    end

    % graph mutant
%             Y = []; E = []; X = [];
%             Y(:,1) = Initial.(M{m}).Y(iG);
%             E(:,1) = Initial.(M{m}).E(iG);
%             X = repmat(tapXlim+2,size(Y));
%             for x = 1:numel(Y)
%                 errorbar(X(x),Y(x),E(x),'MarkerSize',6,'Marker','o',...
%                 'LineStyle','none','LineWidth',1,'Color',color(x+2,1:3),...
%                 'MarkerFaceColor','none',...
%                 'MarkerEdgeColor',color(x+2,1:3));
%             end

    %% put stats star - ref
    if groupN > 1
    p = MWTSet.Graph.Initial.(msr).ANOVA.Prob(1);
    n = MWTSet.Graph.Initial.(msr).Y;
    nd = n(2)-n(1);
    if p < 0.05 
        xx = tapXlim+1;
        yy = rspMax*1.05;
        if nd < 0 
           markertype = 'v';
        elseif nd > 0
           markertype = '^';
        end
        plot(xx,yy,'MarkerSize',6,'Marker',markertype,...
        'LineStyle','none','LineWidth',1,'Color',color(x,1:3),...
        'MarkerFaceColor',color(2,1:3),...
        'MarkerEdgeColor',color(2,1:3));
    end
    end

    %%
%             % put stats star - test
%             a = STAT.(testStrain)(ismember(STAT.msr,M{m}) & ...
%                 ismember(STAT.assays,'Initial'));
%             if strcmp(a,'n.s.') == 0
%                 xx = tapXlim+2;
%                 yy = rspMax*1.05;
%                 if strcmp(a,'increase') == 1
%                    markertype = '^';
%                 elseif strcmp(a,'decrease') == 1
%                    markertype = 'v';
%                 end
%                 plot(xx,yy,'MarkerSize',6,'Marker',markertype,...
%                 'LineStyle','none','LineWidth',1,'Color',color(x,1:3),...
%                 'MarkerFaceColor',color(4,1:3),...
%                 'MarkerEdgeColor',color(4,1:3));
%             end


    %% graph hab dot
%             iG = ~cellfun(@isempty,regexp(Initial.(M{m}).GroupName,['^',strainT]));
    N2i = ~cellfun(@isempty,regexp(Initial.(M{m}).GroupName,'^N2')); % get N2
    % graph N2
    Y = []; E = []; X = [];
    Y(:,1) = HabLevel.(M{m}).Y(N2i);
    E(:,1) = HabLevel.(M{m}).E(N2i);
    X = repmat(tapXlim+1,size(Y));
    for x = 1:numel(Y)
        errorbar(X(x),Y(x),E(x),'MarkerSize',6,'Marker','v',...
        'LineStyle','none','LineWidth',1,'Color',color(x,1:3),...
        'MarkerFaceColor','none',...
        'MarkerEdgeColor',color(x,1:3));
    end

    % graph mutant
%             Y = []; E = []; X = [];
%             Y(:,1) = HabLevel.(M{m}).Y(iG);
%             E(:,1) = HabLevel.(M{m}).E(iG);
%             X = repmat(tapXlim+2,size(Y));
%             for x = 1:numel(Y)
%                 errorbar(X(x),Y(x),E(x),'MarkerSize',6,'Marker','v',...
%                 'LineStyle','none','LineWidth',1,'Color',color(x+2,1:3),...
%                 'MarkerFaceColor','none',...
%                 'MarkerEdgeColor',color(x+2,1:3));
%             end

    % put stats star - ref
    if groupN > 1

    p = MWTSet.Graph.HabLevel.(msr).ANOVA.Prob(1);
    n = MWTSet.Graph.HabLevel.(msr).Y;
    nd = n(2)-n(1);
    if p < 0.05 
        xx = tapXlim+1;
        yy = rspMin-rspMin*.5;
        if nd < 0 
           markertype = 'v';
        elseif nd > 0
           markertype = '^';
        end
        plot(xx,yy,'MarkerSize',6,'Marker',markertype,...
        'LineStyle','none','LineWidth',1,'Color',color(x,1:3),...
        'MarkerFaceColor',color(2,1:3),...
        'MarkerEdgeColor',color(2,1:3));
    end
    end
%             % put stats star - test
%             a = STAT.(testStrain)(ismember(STAT.msr,M{m}) & ...
%                 ismember(STAT.assays,'HabLevel'));
%             if strcmp(a,'n.s.') == 0
%                 xx = tapXlim+2;
%                 yy = rspMin-rspMin*.5;
%                 if strcmp(a,'increase') == 1
%                    markertype = '^';
%                 elseif strcmp(a,'decrease') == 1
%                    markertype = 'v';
%                 end
%                 plot(xx,yy,'MarkerSize',6,'Marker',markertype,...
%                 'LineStyle','none','LineWidth',1,'Color',color(x,1:3),...
%                 'MarkerFaceColor',color(4,1:3),...
%                 'MarkerEdgeColor',color(4,1:3));
%             end


    %% draw divider
    X = repmat(tapN+1,1,2);
    Y = [0 rspMax*1.1];
    line(X,Y,'Color',[0.3 0.3 0.3]);

    %% graph hab rate integral
%             iG = ~cellfun(@isempty,regexp(HabRate.(M{m}).GroupName,['^',strainT]));
    N2i = ~cellfun(@isempty,regexp(HabRate.(M{m}).GroupName,'^N2')); % get N2
%             idata = [find(N2i); find(iG)];
    % graph N2
    Y = []; E = []; 
    % Y will be expressed as average of area under the curve * 3 to
    % fit the graph display
    Y = (HabRate.(M{m}).Y(N2i)./(tapN/1.5))+rspMin; 
    E = (HabRate.(M{m}).E(N2i)./(tapN/1.5));
    X = repmat(tapXlim+4,size(Y));
    for x = 1:numel(Y)
        errorbar(X(x),Y(x),E(x),'MarkerSize',6,'Marker','x',...
        'LineStyle','none','LineWidth',1,'Color',color(x,1:3),...
        'MarkerFaceColor',color(x,1:3),...
        'MarkerEdgeColor',color(x,1:3));
    end
    % graph mutant
%             Y = []; E = []; X = [];
%             Y = HabRate.(M{m}).Y(iG)./(tapN/1.5); 
%             E = HabRate.(M{m}).E(iG)./(tapN/1.5);
%             X = repmat(tapXlim+5,size(Y));
%             for x = 1:numel(Y)
%                 errorbar(X(x),Y(x),E(x),'MarkerSize',6,'Marker','x',...
%                 'LineStyle','none','LineWidth',1,'Color',color(x+2,1:3),...
%                 'MarkerFaceColor',color(x+2,1:3),...
%                 'MarkerEdgeColor',color(x+2,1:3));
%             end
%             

    % put stats star - ref
    if groupN > 1
    p = MWTSet.Graph.HabRate_integral.(msr).ANOVA.Prob(1);
    n = MWTSet.Graph.HabRate_integral.(msr).Y;
    nd = n(2)-n(1);
    if p < 0.05 
        xx = tapXlim+4;
        yy = rspMax*1.05;
        if nd < 0 
           markertype = 'v';
        elseif nd > 0
           markertype = '^';
        end
        plot(xx,yy,'MarkerSize',6,'Marker',markertype,...
        'LineStyle','none','LineWidth',1,'Color',color(x,1:3),...
        'MarkerFaceColor',color(2,1:3),...
        'MarkerEdgeColor',color(2,1:3));
    end
    end

    % put stats star - test
%             a = STAT.(testStrain)(ismember(STAT.msr,M{m}) & ...
%                 ismember(STAT.assays,'HabRate_integral'));
%             if strcmp(a,'n.s.') == 0
%                 xx = tapXlim+5;
%                 yy = rspMax*1.05;
%                 if strcmp(a,'increase') == 1
%                    markertype = '^';
%                 elseif strcmp(a,'decrease') == 1
%                    markertype = 'v';
%                 end
%                 plot(xx,yy,'MarkerSize',6,'Marker',markertype,...
%                 'LineStyle','none','LineWidth',1,'Color',color(x,1:3),...
%                 'MarkerFaceColor',color(4,1:3),...
%                 'MarkerEdgeColor',color(4,1:3));
%             end
%             

    %% draw divider
    X = repmat(tapN+4,1,2);
    Y = [0 rspMax*1.1];
    line(X,Y,'Color',[0.3 0.3 0.3]);

    %% save fig
    savefigepsOnly150([M{m},' ',strain],pSaveA);

end
%     end
% end








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












