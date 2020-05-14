
%% User input
pExpSet = '/Users/connylin/Dropbox/Lab/MWT_Analysis_ByExp/Dose_10sISI';


%% PATHS
% add function path
pFun = ('/Users/connylin/Dropbox/MATLAB/Function_Library_Public');
addpath(pFun);
% add ShaneSpark modules
pProgram = '/Users/connylin/Dropbox/MATLAB/programs_rankinlab/Library/Modules/Dance_ShaneSpark3';
pModules = [pProgram,'/Modules'];
addpath(pModules);
% add module folders
[~,~,~,p] = dircontent(pModules);
for pi = 1:numel(p)
    addpath(p{pi});
end

% create output path
pSave = '/Users/connylin/Dropbox/Lab/Dance Output';
MWTSet.timestamp = generatetimestamp;
pSave = [pSave,'/',MWTSet.timestamp,'_',mfilename];
if isdir(pSave) == 0; mkdir(pSave); end
MWTSet.PATHS.pSaveA = pSave;

% pAnalysis = '/Users/connylin/Dropbox/Lab/MWT_Analysis';
% pData = '/Volumes/ParahippocampalGyrus/MWT/Data';


%% get MWT files from experiment
% get group paths
[~,~,f,p] = dircontent(pExpSet);
% get MWT paths
[~,~,f,p] = cellfun(@dircontent,p,'UniformOutput',0);
pMWTA = celltakeout(p);
pMWTS = pMWTA;

%% IMPORT .TRV % revised 20150412
pMWTf = pMWTA;
A = pMWTA;

% A = MWTSet.MWTInfo.MWTfn_Results;
% pMWTf = MWTSet.MWTInfo.pMWTf_Results;
% MWTfn = MWTSet.MWTInfo.MWTfn_Results;
for m = 1:size(pMWTf,1);
    [~,p] = dircontent(pMWTf{m},'*.trv'); 
    % if there is no .trv
    if isempty(p) == 1
        A{m,2} = {};     
    else       
        % validate trv output format
        pt = p{1};
        fileID = fopen(pt,'r');
        d = textscan(fileID, '%s', 2-1, 'Delimiter', '', 'WhiteSpace', '');
        fclose(fileID);
        % read trv
        if strcmp(d{1}{1}(1),'#') ==1 % if trv file is made by Beethoven
            a = dlmread(pt,' ',5,0); 
        else % if trv file is made by Dance
            a = dlmread(pt,' ',0,0);

        end
        A{m,2} = a(:,[1,3:5,8:10,12:16,19:21,23:27]); % index to none zeros
    end
end
MWTfnImport = A;
MWTSet.Data.Import = MWTfnImport;

% legend
% L = {'time','N?','N_NoResponse','N_Reversed','?','RevDist'    };

%% CHECK TAP CONSISTENCY
[r,c] = cellfun(@size,MWTfnImport(:,2),'UniformOutput',0);
rn = celltakeout(r,'singlenumber');
rfreq = tabulate(rn);
rcommon = rfreq(rfreq(:,2) == max(rfreq(:,2)),1);
str = 'Most common tap number = %d';
display(sprintf(str,rcommon));
rproblem = rn ~= rcommon;
pMWTf = MWTfnImport(:,1);
[~,MWTfn] = cellfun(@fileparts,pMWTf,'UniformOutput',0);


if sum(rproblem)~=0;
    MWTfnP = MWTfn(rproblem); 
    pMWTfP = pMWTf(rproblem);

    str = 'The following MWT did not have the same tap(=%d)';
    display(sprintf(str,rcommon)); 
    disp(MWTfnP);
    % export report
    T = cell2table(MWTfnP,'VariableNames',{'MWT'});
    cd(pSave);
    writetable(T,'MWT_plate_with_missing_taps.csv');
    
    display 'Removing from analysis...'; 
    MWTSet.RawBad = MWTfnImport(rproblem,:);
    MWTfnImport = MWTfnImport(~rproblem,:);
    MWTfnOK = MWTfn(~rproblem);
    pMWTfOK = pMWTf(~rproblem);    

    % reconstruct
%     [MWTSet.MWTfG] = reconstructMWTfG(pMWTf);
    MWTSet.MWTInfo.pMWT = pMWTfOK;
    MWTSet.MWTInfo.MWTfn = MWTfnOK;
    MWTSet.MWTInfo.pMWTBadTap = pMWTf(rproblem); 
    [~,g] = cellfun(@fileparts,...
        cellfun(@fileparts,pMWTfOK,'UniformOutput',0),...
        'UniformOutput',0);
    MWTSet.MWTInfo.GroupName = g;

end

%% MAKING SENSE OF TRV 
% .TRV OUTPUT LEGENDS
% output legends
% time
% # worms already moving backwards (can't score) 
% # worms that didn't reverse in response 
% # worms that did reverse in response 
% mean reversal distance (among those that reversed) 
% standard deviation of reversal distances 
% standard error of the mean 
% minimum 
% 25th percentile 
% median 
% 75th percentile 
% maximum 
% mean duration of reversal (also among those that reversed 
% standard deviation of duration 
% standard error of the mean 
% minimum 
% 25th percentile 
% median 
% 75th percentile 
% maximum

% get data
D = A; 
A = MWTSet.Data.Import(~rproblem,:);
B = struct;
[~,n] = cellfun(@fileparts,A(:,1),'UniformOutput',0);
% B.MWTfn = MWTSet.MWTInfo.MWTfn;
B.MWTfn = n;
B.pMWT = A(:,1);
% pMWTf = MWTSet.MWTInfo.pMWT;
pMWTf = A(:,1);

% indexes of .trv
ind.RevDur = 13;
ind.RevDist = 5;
% calculation

for m = 1:size(pMWTf,1);
    % X = tap time
    % B.X.TapTime(:,m) = A{m,2}(:,1);
    B.X(:,m) = A{m,2}(:,1);   
    % basic caluations
    B.N.NoResponse(:,m) = A{m,2}(:,3);
    B.N.Reversed(:,m) = A{m,2}(:,4);  
    B.N.TotalN(:,m) = B.N.Reversed(:,m)+B.N.NoResponse(:,m);

    %% N
    n = B.N.TotalN(:,m);
    N = B.N.TotalN(:,m);
    N(n < 1) = NaN;

    % Frequency
    B.Y.RevFreq(:,m) = B.N.Reversed(:,m)./N;
    % variance can not be calculated at this point
    B.E.RevFreq(:,m) = NaN(size(B.N.Reversed(:,m))); %  can only be zero
    B.SD.RevFreq(:,m) = NaN(size(B.N.Reversed(:,m)));
    % B.Y.RevFreq(:,m) = B.N.Reversed(:,m)./B.N.TotalN(:,m);
    % B.Y.RevFreqStd(:,m) = B.Y.RevFreq(:,m)/B.Y.RevFreq(1,m);

    % Reversal Duration
    B.Y.RevDur(:,m) = A{m,2}(:,ind.RevDur);
    B.E.RevDur(:,m) = A{m,2}(:,ind.RevDur+1)./B.N.Reversed(:,m);
    B.SD.RevDur(:,m) = A{m,2}(:,ind.RevDur+1)./B.N.Reversed(:,m);

    % Reversal Speed = RevDist/RevDur
    % Distance [disabled]
    RevDist(:,m) = A{m,2}(:,ind.RevDist); 
%     B.SD.RevDist(:,m) = A{m,2}(:,ind.RevDist+1);
%     B.E.RevDist(:,m) = A{m,2}(:,ind.RevDist+1)./B.N.Reversed(:,m);
    % B.Y.RevDistStd(:,m) = B.Y.RevDist(:,m)/B.Y.RevDist(1,m);
    % B.Y.SumRevDist(:,m) = B.Y.RevDist(:,m).*B.N.Reversed(:,m);    
    B.Y.RevSpeed(:,m) = RevDist(:,m)./B.Y.RevDur(:,m); 
    B.E.RevSpeed(:,m) = NaN(size(B.Y.RevSpeed(:,m))); 
    B.SD.RevSpeed(:,m) = NaN(size(B.Y.RevSpeed(:,m)));
end
Raw = B;
MWTSet.Data.Raw = Raw;

%% ORGANIZE OUTPUT BY GROUP
D = MWTSet.Data.Raw;
% GroupName = MWTSet.MWTInfo.GroupName;
pMWT = D.pMWT;
[~,GroupName] = cellfun(@fileparts,...
                cellfun(@fileparts,pMWT,'UniformOutput',0),...
            'UniformOutput',0);
gnameU = unique(GroupName);
% N2 graph first
i = ~cellfun(@isempty,regexp(gnameU(:,1),'^N2'));
% GroupSeq = MWTSet.GraphSetting.GroupSeq;
GroupSeq = gnameU([find(i);find(~i)]);

A = struct;
for g = 1:size(GroupSeq,1)
    B = struct;
    gname = GroupSeq{g};
    i = ismember(GroupName,GroupSeq(g));
    B.MWTplateID = D.MWTfn(i);
    B.time = D.X(:,i);
    B.N_NoResponse = D.N.NoResponse(:,i);
    B.N_Reversed = D.N.Reversed(:,i);
    B.N_TotalN = D.N.TotalN(:,i);
    B.RevFreq_Mean = D.Y.RevFreq(:,i);
    B.RevFreq_SE = D.E.RevFreq(:,i);
    B.RevDur_Mean = D.Y.RevDur(:,i);
    B.RevDur_SE = D.E.RevDur(:,i);
    B.RevSpeed_Mean = D.Y.RevSpeed(:,i);
    B.RevSpeed_SE = D.E.RevSpeed(:,i);
    A.(gname) = B;
end
MWTSet.Data.ByGroupPerPlate = A;

%% MAKE MWT OUTPUT INFO SHEET
[a,MWTfn] = cellfun(@fileparts,pMWTS,'UniformOutput',0);
[a,groupname] = cellfun(@fileparts,a,'UniformOutput',0);
[a,expname] = cellfun(@fileparts,a,'UniformOutput',0);
T = table;
T.MWTfn = MWTfn;
T.groupname = groupname;
T.expname = expname;
cd(pSave);
writetable(T,'expinfosheet.csv');

%% STATISTICS FOR GRAPHS (by group by plate)
D = MWTSet.Data.ByGroupPerPlate;
gnames = fieldnames(D);
B = struct;
B.GroupNames = gnames;
for g = 1:numel(gnames)
    gname = gnames{g};

    D1 = D.(gname);

    plateN = numel(D1.MWTplateID);
    X = (1:size(D1.N_TotalN,1))';

    B.PlateN(g,1) = plateN;

    B.TotalN.X(:,g) = X;
    B.TotalN.Y(:,g) = nanmean(D1.N_TotalN,2);
    B.TotalN.E(:,g) = nanstd(D1.N_TotalN')';

    B.N_Reversed.X(:,g) = X;
    B.N_Reversed.Y(:,g) = nanmean(D1.N_Reversed,2);
    B.N_Reversed.E(:,g) = nanmean(D1.N_Reversed,2);

    B.N_NoResponse.X(:,g) = X;
    B.N_NoResponse.Y(:,g) = nanmean(D1.N_NoResponse,2);
    B.N_NoResponse.E(:,g) = nanstd(D1.N_NoResponse')';

    B.RevFreq.X(:,g) = X;
    B.RevFreq.Y(:,g) = nanmean(D1.RevFreq_Mean,2);
    B.RevFreq.SD(:,g) = nanstd(D1.RevFreq_Mean')';
    B.RevFreq.E(:,g) = B.RevFreq.SD(:,g)./...
        sqrt(repmat(plateN,size(B.RevFreq.SD(:,g),1),1));

    B.RevDur.X(:,g) = X;
    B.RevDur.Y(:,g) = nanmean(D1.RevDur_Mean,2);
    B.RevDur.SD(:,g) = nanstd(D1.RevDur_Mean')';
    B.RevDur.E(:,g) = B.RevDur.SD(:,g)./...
        sqrt(repmat(plateN,size(B.RevDur.SD(:,g),1),1));

    B.RevSpeed.X(:,g) = X;
    B.RevSpeed.Y(:,g) = nanmean(D1.RevSpeed_Mean,2);
    B.RevSpeed.SD(:,g) = nanstd(D1.RevSpeed_Mean')';
    B.RevSpeed.E(:,g) = B.RevSpeed.SD(:,g)./...
        sqrt(repmat(plateN,size(B.RevSpeed.SD(:,g),1),1));

end

MWTSet.Graph.HabCurve = B;

%% GRAPH: HABITUATION CURVES
pSaveA = [MWTSet.PATHS.pSaveA,'/Graph Habituation curves'];
if isdir(pSaveA) == 0; mkdir(pSaveA); end

M = fieldnames(MWTSet.Graph);
M = {'TotalN'
    'N_Reversed'
    'N_NoResponse'
    'RevFreq'
    'RevDur'
    'RevSpeed'};

Graph = MWTSet.Graph.HabCurve;
GroupName = regexprep(Graph.GroupNames,'_',' ');

timestamp = MWTSet.timestamp;
% get N2 
N2i = ~cellfun(@isempty,regexp(GroupName,'^N2'));
% get different strains
% a = regexpcellout(GroupName(~N2i),' ','split');
% strain = unique(a(:,1));

for m = 1:numel(M);% for each measure
    % get coordinates
    X = Graph.(M{m}).X;
    Y = Graph.(M{m}).Y;
    E = Graph.(M{m}).E;

    % create output table legend
    a = cellstr([repmat('t',size(Y,1),1) num2str([1:size(Y,1)]')]);
    a = regexprep(a,' ','');
    b = cellstr([repmat('SE',size(Y,1),1) num2str([1:size(Y,1)]')]);
    b = regexprep(b,' ','');
    vnames = [{'groupname'}; a;b];   
    % make table 
    d = [GroupName num2cell(Y') num2cell(E')];
    T = cell2table(d,'VariableNames',vnames);
    cd(pSaveA);
    writetable(T,[M{m},'.csv'],'Delimiter',',');
        
    % Create figure
    figure1 = figure('Color',[1 1 1],'Visible','off');
    axes1 = axes('Parent',figure1,'FontSize',18);
    box(axes1,'off');
    xlim(axes1,[0 size(Y,1)+1]);
    ylim(axes1,[min(min(Y-E))*.9 max(max(Y+E))*1.1]);
    hold(axes1,'all');
    errorbar1 = errorbar(X,Y,E,'LineWidth',1);
    color = [0 0 0; 1 0 0; [0.5 0.5 0.5]; [0.04 0.52 0.78]]; 
    n = min([size(color,1) size(Y,2)]);
    for x = 1:n
        set(errorbar1(x),'Color',color(x,1:3))
    end

    for x = 1:size(Y,2)
        set(errorbar1(x),'DisplayName',GroupName{x})
    end
    titlestr = timestamp;
    title(titlestr,'FontSize',12);
    xlabel('Tap','FontSize',18);
    ylabel(regexprep(M{m},'_',' '),'FontSize',18);
    legend1 = legend(axes1,'show');
    set(legend1,'EdgeColor',[1 1 1],...
        'Location','NorthEastOutside',...
        'YColor',[1 1 1],'XColor',[1 1 1],'FontSize',12);

    % text strings
    a = Graph.PlateN;
    b = '';
    for x = 1:numel(a)
        b = [b,num2str(a(x)),', '];
    end
    str1 = ['N(plate) = ',b(:,1:end-2)];
    textboxstr = [str1];
    annotation(figure1,'textbox',...
        [0.70 0.015 0.256 0.05],...
        'String',textboxstr,...
        'FitBoxToText','off',...
        'EdgeColor',[1 1 1]);
    % save fig
    savefigepsOnly(M{m},pSaveA);
    
end

%% GRAPH: HABITUATION CURVES BY STRAINS
pSaveA = [MWTSet.PATHS.pSaveA,'/Graph Habituation curves by strain'];
if isdir(pSaveA) == 0; mkdir(pSaveA); end

M = fieldnames(MWTSet.Graph);
M = {'TotalN'
    'N_Reversed'
    'N_NoResponse'
    'RevFreq'
    'RevDur'
    'RevSpeed'};

Graph = MWTSet.Graph.HabCurve;
GroupName = regexprep(Graph.GroupNames,'_',' ');

timestamp = MWTSet.timestamp;
% get N2 
N2i = ~cellfun(@isempty,regexp(GroupName,'^N2'));
if sum(N2i) ~= numel(GroupName)
    % get different strains
    a = regexpcellout(GroupName(~N2i),' ','split');
    strain = unique(a(:,1));

    for m = 1:numel(M);% for each measure
        % get coordinates
        X = Graph.(M{m}).X;
        Y = Graph.(M{m}).Y;
        E = Graph.(M{m}).E;

        % make figures by group and N2
        X1 = X; Y1 = Y;  E1 = E;
        for straini = 1:numel(strain)
            strainT = strain{straini};
            pSaveA2 = [pSaveA,'/',strainT];
            if isdir(pSaveA2) == 0; mkdir(pSaveA2); end
            iG = ~cellfun(@isempty,regexp(GroupName,['^',strainT]));
            gname = GroupName([find(N2i);find(iG)]);
            X = repmat(X1(:,1),1,sum(N2i)+sum(iG)); 
            Y = Y1(:,[find(N2i);find(iG)]); 
            E = E1(:,[find(N2i);find(iG)]);

            % create figure
            figure1 = figure('Color',[1 1 1],'Visible','off');
            axes1 = axes('Parent',figure1,'FontSize',24);
            box(axes1,'off');
            xlim(axes1,[0 size(Y,1)+1]);
            ylim(axes1,[min(min(Y-E))*.9 max(max(Y+E))*1.1]);
            hold(axes1,'all');
            errorbar1 = errorbar(X,Y,E,'LineWidth',1);
            color = [0 0 0; 1 0 0; [0.5 0.5 0.5]; [0.04 0.52 0.78]]; 
            n = min([size(color,1) size(Y,2)]);
            for x = 1:n
                set(errorbar1(x),'Color',color(x,1:3))
            end

            for x = 1:size(Y,2)
                set(errorbar1(x),'DisplayName',gname{x})
            end
            titlestr = timestamp;
            title(titlestr,'FontSize',12);
            xlabel('Tap','FontSize',18);
            ylabel(regexprep(M{m},'_',' '),'FontSize',18);
            legend1 = legend(axes1,'show');
            set(legend1,'EdgeColor',[1 1 1],...
                'Location','NorthEastOutside',...
                'YColor',[1 1 1],'XColor',[1 1 1],'FontSize',12);

            % text strings
            a = Graph.PlateN(iG);
            b = '';
            for x = 1:numel(a)
                b = [b,num2str(a(x)),', '];
            end
            str1 = ['N(plate) = ',b(:,1:end-2)];
            textboxstr = [str1];
            annotation(figure1,'textbox',...
                [0.70 0.015 0.256 0.05],...
                'String',textboxstr,...
                'FitBoxToText','off',...
                'EdgeColor',[1 1 1]);
            % save fig
            savefigepsOnly150([M{m},' ',strain{straini}],pSaveA2);
    %         savefigeps([M{m},' ',strain{straini}],pSaveA2);

        end
    end
end

%% EXCEL OUTPUT -GROUP
Data = MWTSet.Graph.HabCurve;
GroupNames = MWTSet.Graph.HabCurve.GroupNames;
MsrList = {'RevFreq','RevDur','RevSpeed'};
for Msri = 1:numel(MsrList)
    D = Data.(MsrList{Msri});
    d = [GroupNames num2cell(D.Y') num2cell(D.E')];

    a = cellstr([repmat('t',size(D.Y,1),1) num2str([1:size(D.Y,1)]')]);
    a = regexprep(a,' ','');
    b = cellstr([repmat('SE',size(D.Y,1),1) num2str([1:size(D.Y,1)]')]);
    b = regexprep(b,' ','');
    vnames = [{'groupname'}; a;b];
    
    T = cell2table(d,'VariableNames',vnames);
    cd(pSaveA);
    writetable(T,[MsrList{Msri},'.dat'],'Delimiter','\t');
end

%% CALCULATION: LAST 3 TAPS 
pSaveA = [pSave,'/','Stats HabLevel'];
if isdir(pSaveA) == 0; mkdir(pSaveA); end
M = fieldnames(MWTSet.Graph);
M = {'RevFreq', 'RevDur', 'RevSpeed'};

% calculate data - hab level
D = MWTSet.Data.ByGroupPerPlate;
gnames = fieldnames(D);
A = struct;
for m = 1:numel(M);% for each measure
    B = []; G = {};
    for g = 1:numel(gnames)
        d = D.(gnames{g}).([M{m},'_Mean'])(end-2:end,:);
        % calculate last 3 taps average per plate
        d1 = nanmean(d);
        % get data
        B = [B;d1'];
        n = numel(d1);
        G = [G; repmat(gnames(g), n,1)];
    end
    % stats
    [m2, n2, se2,gnames2] = grpstats(B,G,{'mean','numel','sem','gname'});
    A.(M{m}).GroupName = gnames2;
    A.(M{m}).N = n2;
    A.(M{m}).Y = m2;
    A.(M{m}).E = se2; 

    if numel(unique(G)) > 1
        % anova
        [p,t,stats] = anova1(B,G,'off');
        [c,m1,h,gnames] = multcompare(stats,'ctype','bonferroni','display','off');
        A.(M{m}).ANOVA = t;
        % export anova
        vname = regexprep(t(1,:),'Prob>F','P_value');
        T = cell2table(t(2:end,:),'VariableNames',vname);
        cd(pSaveA); 
        writetable(T,[M{m},' HabLevel ANOVA.csv'],'Delimiter',',');
        
        % posthoc
        i = ((c(:,2) >0 & c(:,4) >0) + (c(:,2) <0 & c(:,4) <0)) >0;
        a = [gnames(c(:,1)), gnames(c(:,2))];
        a(i,3) = {'< 0.05'};
        a(~i,3) = {'n.s.'};
        A.(M{m}).posthoc = a;
        % export posthoc
        T = cell2table(a);
        cd(pSaveA); 
        writetable(T,[M{m},' HabLevel posthoc bonferroni.csv'],'Delimiter',',');
        
    end
end
MWTSet.Graph.HabLevel = A;

%% CALCULATION: INITIAL
pSaveA = [pSave,'/','Stats Initial'];
if isdir(pSaveA) == 0; mkdir(pSaveA); end
D = MWTSet.Data.ByGroupPerPlate;
gnames = fieldnames(D);
A = struct;
for m = 1:numel(M);% for each measure
    B = []; G = {};
    for g = 1:numel(gnames)
        d1 = D.(gnames{g}).([M{m},'_Mean'])(1,:);
        % get data
        B = [B;d1'];
        n = numel(d1);
        G = [G; repmat(gnames(g), n,1)];

    end

    % stats
    [m2, n2, se2,gnames2] = grpstats(B,G,{'mean','numel','sem','gname'});
    A.(M{m}).GroupName = gnames2;
    A.(M{m}).N = n2;
    A.(M{m}).Y = m2;
    A.(M{m}).E = se2; 

    if numel(unique(G)) > 1
        % anova
        [p,t,stats] = anova1(B,G,'off');
        [c,m1,h,gnames] = multcompare(stats,'ctype','bonferroni','display','off');
        A.(M{m}).ANOVA = t;
        % export anova
        vname = regexprep(t(1,:),'Prob>F','P_value');
        T = cell2table(t(2:end,:),'VariableNames',vname);
        cd(pSaveA); 
        writetable(T,[M{m},' Initial ANOVA.csv'],'Delimiter',',');
        
        % posthoc
        i = ((c(:,2) >0 & c(:,4) >0) + (c(:,2) <0 & c(:,4) <0)) >0;
        a = [gnames(c(:,1)), gnames(c(:,2))];
        a(i,3) = {'< 0.05'};
        a(~i,3) = {'not significant'};
        A.(M{m}).posthoc = a;
        % export posthoc
        T = cell2table(a);
        cd(pSaveA); 
        writetable(T,[M{m},' Initial posthoc bonferroni.csv'],'Delimiter',',');
        
    end
end
MWTSet.Graph.Initial = A;

%% GRAPH: INITIAL HAB LEVEL BAR
pSaveA = [pSave,'/','Graph Initial HabLevel bar'];
if isdir(pSaveA) == 0; mkdir(pSaveA); end

gcategory = {'Initial','HabLevel'};
% load function module
addpath(sprintf('%s/xticklabel_rotate',pFun));

G = []; 
for m = 1:numel(M); % for each measure
    GroupName = {}; Y = []; E = []; X = [];
    for c = 1:numel(gcategory)
        GroupName(:,c) = MWTSet.Graph.(gcategory{c}).(M{m}).GroupName;
        X(:,c) = (1:numel(MWTSet.Graph.(gcategory{c}).(M{m}).Y))';
        Y(:,c) = MWTSet.Graph.(gcategory{c}).(M{m}).Y;
        E(:,c) = MWTSet.Graph.(gcategory{c}).(M{m}).E;
    end
    
    % create output table legend
    vnames = {'groupname', 'Initial','HabLevel','Initial_SE','HabLevel_SE'};   
    % make table 
    d = [GroupName(:,1) num2cell(Y) num2cell(E)];
    T = cell2table(d,'VariableNames',vnames);
    cd(pSaveA);
    writetable(T,[M{m},'Initial HabLevel.csv'],'Delimiter',',');    
  
    % GRAPHING
    titlename = MWTSet.timestamp;
    figname  = sprintf('%s Initial HabLevel bar',M{m});
%     figure1 = figure('Color',[1 1 1],'Visible','off');
    figure1 = figure('Color',[1 1 1],'Visible','off');
    axes1 = axes('Parent',figure1,...
        'XTickLabel',regexprep(GroupName(:,1),'_',' '),...
        'XTick',1:size(Y,1),...
        'FontSize',12);
    xlim(axes1,[0.5 size(Y,1)+0.5]);
    hold(axes1,'all');

    % create bar
    bar1 = bar(X,Y,'BarWidth',1,'Parent',axes1); 
    colorSet = [[0.043 0.52 0.78];[0.85 0.16 0]];
    for c = 1:size(colorSet,1)
        set(bar1(c),'FaceColor',colorSet(c,:),'DisplayName',gcategory{c});
    end

    % create errorbar
    errorbar1 = errorbar([X(:,1)-0.145 X(:,2)+0.145],Y,E);
    for c = 1:size(colorSet,1)
        set(errorbar1(c),...
            'LineStyle','none','Color',[0 0 0],'DisplayName',' ');
    end

    title(titlename);
    ylabel(M{m},'FontSize',16);
    legend1 = legend(axes1,'show');
    set(legend1,'EdgeColor',[1 1 1],...
        'Location','EastOutside','YColor',[1 1 1],'XColor',[1 1 1]);
    xticklabel_rotate; % rotate X labels
%     savefigepsOnly150(figname,pSaveA);   
    savefigepsOnly(figname,pSaveA);    
end

%% GRAPH: INITAL VS HAB LEVEL DOT GRAPH
% if size(MWTSet.GraphSetting.GroupSeq,1) > 1

    pSaveA = [MWTSet.PATHS.pSaveA,'/','Graph Initial HabLevel dot'];
    if isdir(pSaveA) == 0; mkdir(pSaveA); end
    M = {'RevFreq', 'RevDur', 'RevSpeed'};
    for m = 1:numel(M)
        Y = []; E = []; X = [];
        gname = regexprep(MWTSet.Graph.Initial.(M{m}).GroupName,'_',' ');
        Y(:,1) = MWTSet.Graph.Initial.(M{m}).Y;
        E(:,1) = MWTSet.Graph.Initial.(M{m}).E;
        Y(:,2) = MWTSet.Graph.HabLevel.(M{m}).Y;
        E(:,2) = MWTSet.Graph.HabLevel.(M{m}).E;
        X(:,1:size(Y,2)) = repmat((1:size(Y,1))',1,size(Y,2));
    
        % make graph
        figure1 = figure('Color',[1 1 1],'Visible','off');
        axes1 = axes('Parent',figure1,...
            'XTickLabel',gname,...
            'XTick',1:size(Y,1),...
            'FontSize',14);
        hold(axes1,'all');
        errorbar1 = errorbar(X,Y,E,...
            'MarkerSize',6,'Marker','o',...
            'LineStyle','none',...
            'LineWidth',1);
        set(errorbar1(1),'MarkerFaceColor',[0 0 0],...
            'MarkerEdgeColor',[0 0 0],...
            'DisplayName','Initial',...
            'Color',[0 0 0]);
        set(errorbar1(2),'MarkerFaceColor',[1 0 0],...
            'MarkerEdgeColor',[1 0 0],...
            'DisplayName','HabLevel',...
            'Color',[1 0 0]);
        ylabel(M{m},'FontSize',18);
        legend1 = legend(axes1,'show');
        set(legend1,'Location','EastOutside','EdgeColor',[1 1 1],'YColor',[1 1 1],'XColor',[1 1 1]);
%         savefigepsOnly150(['Initial vs hab level',M{m}],pSaveA);
        xticklabel_rotate; % rotate X labels
        savefigepsOnly([M{m},' Initial vs hab level'],pSaveA);    

    end
% end

%% GRAPH: HAB STRENGTH (INITIAL - HAB LEVEL)
pSaveA = [MWTSet.PATHS.pSaveA,'/','Graph Initial minus HabLevel'];
if isdir(pSaveA) == 0; mkdir(pSaveA); end
D = MWTSet.Data.ByGroupPerPlate;
m = 1;
gnames = fieldnames(D);
for m = 1:numel(M)
    B = [];
    G = [];
    for g = 1:numel(gnames)
        % initial
        % hab level
        d = D.(gnames{g}).([M{m},'_Mean'])(end-2:end,:);
        hablevel = nanmean(d);
        initial = D.(gnames{g}).([M{m},'_Mean'])(1,:);
        d1 = initial - hablevel;
        % calculate last 3 taps average per plate 
        % get data
        B = [B;d1'];
        n = numel(d1);
        G = [G; repmat(gnames(g), n,1)];
    end
    % stats
    [m2, n2, se2,gnames2] = grpstats(B,G,{'mean','numel','sem','gname'});
    A.(M{m}).GroupName = gnames2;
    A.(M{m}).N = n2;
    A.(M{m}).Y = m2;
    A.(M{m}).E = se2; 

    if numel(unique(G)) > 1
        % anova
        [p,t,stats] = anova1(B,G,'off');
        [c,m1,h,gnames] = multcompare(stats,'ctype','bonferroni','display','off');
        A.(M{m}).ANOVA = t;
        % export anova
        vname = regexprep(t(1,:),'Prob>F','P_value');
        T = cell2table(t(2:end,:),'VariableNames',vname);
        cd(pSaveA); 
        writetable(T,[M{m},' Hab Strength ANOVA.csv'],'Delimiter',',');

        i = ((c(:,2) >0 & c(:,4) >0) + (c(:,2) <0 & c(:,4) <0)) >0;
        a = [gnames(c(:,1)), gnames(c(:,2))];
        a(i,3) = {'< 0.05'};
        a(~i,3) = {'not significant'};
        A.(M{m}).posthoc = a;
        % export posthoc
        T = cell2table(a);
        cd(pSaveA); 
        writetable(T,[M{m},' Hab Strength posthoc bonferroni.csv'],'Delimiter',',');
    end
end
MWTSet.Graph.HabStrength = A;

% graphing
Y = []; X = []; E = [];
for m = 1:numel(M)
    gname = regexprep(MWTSet.Graph.HabStrength.(M{m}).GroupName,'_',' ');
    Y = MWTSet.Graph.HabStrength.(M{m}).Y;
    E = MWTSet.Graph.HabStrength.(M{m}).E;
    X(:,1:size(Y,2)) = repmat((1:size(Y,1))',1,size(Y,2));
    
    
    % graph
    figure1 = figure('Color',[1 1 1],'Visible','off');
    axes1 = axes('Parent',figure1,...
        'XTickLabel',gname,...
        'XTick',1:size(Y,1),...
        'FontSize',14);
    hold(axes1,'all');
    errorbar1 = errorbar(X,Y,E,'MarkerSize',10,'Marker','o',...
        'LineStyle','none',...
        'LineWidth',1);
    set(errorbar1(1),'MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[0 0 0],...
        'DisplayName','Initial',...
        'Color',[0 0 0]);

    ylabel([M{m}, ' (initial - hab level)'],'FontSize',18);
%     legend1 = legend(axes1,'show');
%     set(legend1,'Location','EastOutside','EdgeColor',[1 1 1],'YColor',[1 1 1],'XColor',[1 1 1]);
%     savefigepsOnly150(['Habituation stength ',M{m}],pSaveA);
    xticklabel_rotate; % rotate X labels
    savefigepsOnly([M{m},' Hab strength'],pSaveA);  
end

%% GRAPH: HAB STRENGTH % (percent decrease Hab level - INITIAL)
pSaveA = [MWTSet.PATHS.pSaveA,'/','Graph Initial vs HabLevel'];
if isdir(pSaveA) == 0; mkdir(pSaveA); end
D = MWTSet.Data.ByGroupPerPlate;
m = 1;
gnames = fieldnames(D);
for m = 1:numel(M)
    B = [];
    G = [];
    for g = 1:numel(gnames)
        % hab level
        d = D.(gnames{g}).([M{m},'_Mean'])(end-2:end,:);
        hablevel = nanmean(d);
        % initial
        initial = D.(gnames{g}).([M{m},'_Mean'])(1,:);
        d1 = ((hablevel - initial)./initial)*100;
        % calculate last 3 taps average per plate 
        % get data
        B = [B;d1'];
        n = numel(d1);
        G = [G; repmat(gnames(g), n,1)];
        B(B == Inf) = nan;
    end
        % stats
    [m2, n2, se2,gnames2] = grpstats(B,G,{'mean','numel','sem','gname'});
    A.(M{m}).GroupName = gnames2;
    A.(M{m}).N = n2;
    A.(M{m}).Y = m2;
    A.(M{m}).E = se2; 

    if numel(unique(G)) > 1
        % anova
        [p,t,stats] = anova1(B,G,'off');
        [c,m1,h,gnames] = multcompare(stats,'ctype','bonferroni','display','off');
        A.(M{m}).ANOVA = t;

        i = ((c(:,2) >0 & c(:,4) >0) + (c(:,2) <0 & c(:,4) <0)) >0;
        a = [gnames(c(:,1)), gnames(c(:,2))];
        a(i,3) = {'< 0.05'};
        a(~i,3) = {'not significant'};
        A.(M{m}).posthoc = a;
    end
end
MWTSet.Graph.HabStrength_percent = A;

% graphing
Y = []; X = []; E = [];

errorflag = 0;

for m = 1:numel(M)
    gnnow = MWTSet.Graph.HabStrength_percent.(M{m}).GroupName;
    if m > 1 && numel(gnnow) ~= size(Y,1)
        warning('some groups has missing %s data',M{m});
        errorflag = 1;
    else
        Y(:,m) = MWTSet.Graph.HabStrength_percent.(M{m}).Y;
        E(:,m) = MWTSet.Graph.HabStrength_percent.(M{m}).E;
        X(:,1:size(Y,2)) = repmat((1:size(Y,1))',1,size(Y,2));
    end
end    

gname = gnnow;

if errorflag == 0;
    figure1 = figure('Color',[1 1 1],'Visible','off');
    axes1 = axes('Parent',figure1,...
        'XTickLabel',regexprep(gname,'_',' '),...
        'XTick',1:size(Y,1),...
        'FontSize',12);
    hold(axes1,'all');
    errorbar1 = errorbar(X,Y,E,'MarkerSize',6,'Marker','o',...
        'LineStyle','none',...
        'LineWidth',1);
    colorset = [0 0 0; 1 0 0; 0 0 1];
    if size(colorset,1) > size(X,2)
        xN = size(X,2);
    else
        xN = size(colorset,1);
    end
    for x = 1:xN
        set(errorbar1(x),...
            'MarkerFaceColor',colorset(x,:),...
            'MarkerEdgeColor',colorset(x,:),...
            'DisplayName',M{x},...
            'Color',colorset(x,:));
    end

    ylabel('% diff (hab level - initial)','FontSize',16);
    legend1 = legend(axes1,'show');
    set(legend1,'Location','NorthOutside','EdgeColor',[1 1 1],'YColor',[1 1 1],'XColor',[1 1 1]);
%     savefigepsOnly150(['Habituation stength percent difference',M{m}],pSaveA);
    xticklabel_rotate; % rotate X labels
    savefigepsOnly('Hab stength percent difference',pSaveA);  
else
    warning('no [Habituation stength percent difference] graph created');
end



%% EXCEL OUTPUT
pSaveA = [MWTSet.PATHS.pSaveA,'/','Tables'];
if isdir(pSaveA) == 0; mkdir(pSaveA); end
gcategory = {'Initial','HabLevel'};
groupname = {};
rtype = {};
ctype = {};
N = []; Y = []; E = [];
for c = 1:numel(gcategory)
    for m = 1:numel(M)
        D = MWTSet.Graph.(gcategory{c}).(M{m});
        groupname = [groupname; D.GroupName];
        rtype = [rtype; repmat(M(m), numel(D.GroupName), 1)];
        ctype = [ctype; repmat(gcategory(c), numel(D.GroupName), 1)];
        N = [N; D.N];
        Y = [Y; D.Y];
        E = [E; D.E];       
    end
end
T = table;
T.group = groupname;
T.response_type = rtype;
T.msr = ctype;
T.N_plates = N;
T.Mean = Y;
T.SE = E;

cd(pSaveA);
writetable(T,'initial hab level summary.txt','Delimiter','\t');

% export anova 
for m = 1:numel(M)
    for c = 1:numel(gcategory)
        nn = fieldnames(MWTSet.Graph.(gcategory{c}).(M{m}));
        if strcmp(nn,'AVOVA') == 1
            D = MWTSet.Graph.(gcategory{c}).(M{m}).ANOVA;
            D(1,6) = {'p'};
            T = cell2table(D(2:end,:),'VariableNames',D(1,:));
            tname = sprintf('ANOVA %s %s.txt',gcategory{c},M{m});
            cd(pSaveA); 
            writetable(T, tname,'Delimiter','\t');
        end
    end   
end
% export posthoc
for m = 1:numel(M)
    for c = 1:numel(gcategory)
        nn = fieldnames(MWTSet.Graph.(gcategory{c}).(M{m}));
        if strcmp(nn,'posthoc') == 1
            D = MWTSet.Graph.(gcategory{c}).(M{m}).posthoc;
            T = cell2table(D,'VariableNames',{'group1','group2','p'});
            tname = sprintf('ANOVA %s %s posthoc bonferroni.txt',...
                gcategory{c},M{m});
            cd(pSaveA); 
            writetable(T, tname,'Delimiter','\t');
        end
    end   
end

%% SAVE MAT FILES
cd(pSave);
save('matlab.mat')


%% Report done
fprintf('\n\n***DONE***\n\n');


