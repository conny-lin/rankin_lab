function Stats_HabCurveRMANOVA(Data,MWTDB,pSave,varargin)

%% default
alpha = 0.05;
pvlimit = 0.001;
compName = 'groupname';
prefix = '';
suffix = '';
msrlist = {'RevFreq','RevSpeed','RevDur'};
    
%%
vararginProcessor


%% create input table
[DataS,tptable] = DataPrep_rANOVA(Data);

DataS.groupname = MWTDB.groupname(DataS.mwtid);
DataS.strain = MWTDB.strain(DataS.mwtid);
DataS.dose = num2cellstr(parseGname_TestmM(DataS.groupname));



% tpu = unique(Data.tap);
% measNames = num2cellstr(tpu);
% tpnames = strjoinrows([cellfunexpr(measNames,'tap'),measNames],'');
% 
% 
% mwtidu = unique(Data.mwtid);
% A = nan(numel(mwtidu)*numel(msrlist),numel(tpu));
% Leg = cell(numel(mwtidu)*numel(msrlist),2);
% row1 = 1;
% for mwti = 1:numel(mwtidu)
%     for msri = 1:numel(msrlist)
%         d = table2array(Data(Data.mwtid==mwtidu(mwti),msrlist(msri)));
%         col = table2array(Data(Data.mwtid==mwtidu(mwti),{'tap'}));
%         A(row1,col) = d;
%         Leg{row1,1} = msrlist{msri};
%         Leg{row1,2} = mwtidu(mwti);
%         row1 = row1+1;
%     end
% end
% A = array2table(A,'VariableNames',tpnames);
% Leg = cell2table(Leg,'VariableNames',{'msr','mwtid'});
% Leg.groupname = MWTDB.groupname(Leg.mwtid);
% Leg.strain = MWTDB.strain(Leg.mwtid);
% Leg.dose = num2cellstr(parseGname_TestmM(Leg.groupname));
% DataS = [Leg A];
% 
% tptable = table(tpu,'VariableNames',{'tap'});


%% create rmanova options

multiStrain = numel(unique(DataS.strain)) >1;
multiDose = numel(unique(DataS.dose)) > 1;

gpairs = strjoinrows(pairwisecomp_getpairs(unique(DataS.(compName))),' x ');

%%

cd(pSave); close all;
filesavename = sprintf('%sRMANOVA%s.txt',prefix,suffix);
fid = fopen(filesavename,'w');
for msri = 1:numel(msrlist)
    msr= msrlist{msri};
    if msri==1
        fprintf(fid,'\n----- %s -----\n',msr);
    else
        fprintf(fid,'\n\n----- %s -----\n',msr);
    end
    %% get data
    i = ismember(DataS.msr,msr);
    D = DataS(i,:);
    
    % decide which to run
    if multiStrain && multiDose; factorName = 'strain*dose';
    elseif multiStrain && ~multiDose; factorName = 'strain';
    elseif ~multiStrain && multiDose; factorName = 'dose';
    end
    rmTerms = sprintf('tap1-tap30~%s',factorName);
    rm = fitrm(D,rmTerms,'WithinDesign',tptable);
    fprintf(fid,'RMANOVA(tap:%s):\n%s\n',factorName,anovan_textresult(ranova(rm), 0, 'pvlimit',pvlimit));

    % run pairwise
    rm = fitrm(D,'tap1-tap30~groupname','WithinDesign',tptable);
    t = multcompare(rm,'groupname');
    
    
    fprintf(fid,'\nPosthoc(Tukey) HabCurve by %s:\n',compName);
    if isempty(t); fprintf(fid,'All comparison = n.s.\n');
    else fprintf(fid,'%s\n',multcompare_text2016b(t,'pvlimit',pvlimit,'alpha',alpha));
    end
    
    % comparison by taps
    t = multcompare(rm,compName,'By','tap');
    %  keep only unique comparisons
    a = strjoinrows([t.([compName,'_1']) t.([compName,'_2'])],' x ');
    t(~ismember(a,gpairs),:) =[]; 
    % record
    fprintf(fid,'\nPosthoc(Tukey)tap by %s:\n',compName); 
    if isempty(t); fprintf(fid,'All comparison = n.s.\n');
    else fprintf(fid,'%s\n',multcompare_text2016b(t,'pvlimit',pvlimit,'alpha',alpha));    
    end
end

fclose(fid);















