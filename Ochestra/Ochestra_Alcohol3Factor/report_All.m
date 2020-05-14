function report_All(pData_Sen,pData_TAR,pData_SS,pSave,MWTDB,strain)
%% STATS REPORT
%% curve sensitivity +++++++++++++++
curveData = rmANOVACurveReport;
curveData.datapath = fullfile(pData_Sen,'curve MANOVA.txt');
%% calculate percentage difference 
pC = fullfile(pData_Sen,'data.mat');
DC = load(pC,'DataMeta','MWTDB');
CS = CurveStats;
CS.mwtid = DC.DataMeta.mwtid;
CS.curve = DC.DataMeta.curve;
CS.MWTDB = DC.MWTDB;
txt = reportCurveSensitivity(curveData,CS);
%% ---------------------------------

% TAR ++++++++++
TAR = rmANOVATARReport;
TAR.datapath = fullfile(pData_TAR, 'AccProb RMANOVA.txt');
txt2 = reportTAR(TAR);
% --------------

% TWR  ++++++++++
TWR = rmANOVATWR;
TWR.datapath = fullfile(pData_SS, 'RMANOVA.txt');
txt3 = reportTWR(TWR);
% --------------

% N ++++++++++
g = curveData.groupnames';
n = curveData.n';
n = sortN2first(g,n);
n = strjoin(num2cellstr(curveData.n),',');
g = sortN2first(g,g);
g = strjoin(g,', ');
nt1 = sprintf('Curve, N(worms, %s) = %s',g,n);
t = tabulate(MWTDB.groupname);
gname = sortN2first(t(:,1),t(:,1));
gname = regexprep(gname,'_',' ');
gname = strjoin(gname,', ');
n = sortN2first(t(:,1),t(:,2));
n = strjoin(cellfun(@num2cellstr,n),', ');
ntxt = sprintf('%s\nResponses, N(plates, %s) = %s',nt1,gname,n);
% --------------

txtFinal = sprintf('%s\n%s\n%s\n\n',txt,txt2,txt3,ntxt);


% EXPORT REPORT +++++++++
p2 = fullfile(pSave,sprintf('%s stats writeup.txt',strain));
fid = fopen(p2,'w');
fprintf(fid,'%s',txtFinal);
fclose(fid);
%-----------

% save objects ++++
p3 = fullfile(pSave,sprintf('%s.mat',strain));
save(p3,'curveData','TAR','TWR','MWTDB');
% ---------------