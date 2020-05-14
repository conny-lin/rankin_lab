function Stats_HabCurveRMANOVA_v2(Data,MWTDB,pSave,varargin)

%% default
alpha = 0.05;
pvlimit = 0.001;
compName = 'groupname';
prefix = '';
suffix = '';
space = true;
msrlist = {'RevFreq','RevSpeed','RevDur'};
    
%% varargin
vararginProcessor


%% create input table
[DataS,tptable] = DataPrep_rANOVA(Data);
DataS.groupname = MWTDB.groupname(DataS.mwtid);
DataS.strain = MWTDB.strain(DataS.mwtid);
DataS.dose = num2cellstr(parseGname_TestmM(DataS.groupname));


%% create rmanova options
multiStrain = numel(unique(DataS.strain)) >1;
multiDose = numel(unique(DataS.dose)) > 1;
gpairs = strjoinrows(pairwisecomp_getpairs(unique(DataS.(compName))),' x ');

%% generate output
filesavename = fullfile(pSave,sprintf('%sRMANOVA%s.txt',prefix,suffix));
fid = fopen(filesavename,'w');
for msri = 1:numel(msrlist)
    msr= msrlist{msri}; % get msr name
    
    % write msr title
    if msri==1
        fprintf(fid,'\n----- %s -----\n',msr);
    else
        fprintf(fid,'\n\n----- %s -----\n',msr);
    end
    
    %% descriptive stats
    a = tabulate(DataS.groupname(ismember(DataS.msr,msrlist{1}),:));
    gnames = a(:,1);
    if space; s = ', '; else; s = ','; end
    N = strjoin(num2cellstr(cell2mat(a(:,2)))',s);
    gstr = strjoin(gnames',s);
    if space; s = ' '; else; s = ''; end
    Nstr = sprintf('N%s=%s%s',s, s, N);
    fprintf(fid,'%s\n%s\n',gstr,Nstr);
    %%
    % get data
    i = ismember(DataS.msr,msr);
    D = DataS(i,:);
    
    % decide which to run
    if multiStrain && multiDose
           factorName = 'strain*dose';
    elseif multiStrain && ~multiDose
        factorName = 'strain';
    elseif ~multiStrain && multiDose
        factorName = 'dose';
    end
    
    % rmanova
    rmTerms = sprintf('tap1-tap30~%s',factorName); % get anova terms
    rm = fitrm(D,rmTerms,'WithinDesign',tptable); % run rmanova
   
    % print results
    ranovatable = ranova(rm);
    txt = anovan_textresult(ranovatable, 0, 'pvlimit',pvlimit,'Fdigit',2,'pvspace',1); % get anova text result
    fprintf(fid,'RMANOVA(tap*%s):\n%s\n',factorName,txt); % write results

    % run pairwise
    rm = fitrm(D,'tap1-tap30~groupname','WithinDesign',tptable);
    t = multcompare(rm,'groupname');
    
    fprintf(fid,'\nPosthoc(Tukey) HabCurve by %s:\n',compName);
    if isempty(t)
        fprintf(fid,'All comparison = n.s.\n');
    else
        fprintf(fid,'%s\n',multcompare_text2016b(t,'pvlimit',pvlimit,'alpha',alpha));
    end
    
    % comparison by taps
    t = multcompare(rm,compName,'By','tap');
    
    %  keep only unique comparisons
    a = strjoinrows([t.([compName,'_1']) t.([compName,'_2'])],' x ');
    t(~ismember(a,gpairs),:) =[]; 
    
    % record
    fprintf(fid,'\nPosthoc(Tukey)tap by %s:\n',compName); 
    if isempty(t)
        fprintf(fid,'All comparison = n.s.\n');
    else
        fprintf(fid,'%s\n',multcompare_text2016b(t,'pvlimit',pvlimit,'alpha',alpha));
    end
end

fclose(fid);















