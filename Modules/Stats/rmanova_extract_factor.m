function [TM, DVnames,IVnames, dmTable, dmFactorNum, T2] = rmanova_extract_factor(T,varargin)

%% default
rmfactorName = 't';
factor1Name = 'gname';
groupnamelist= {''};
vararginProcessor;



%% prep for rmanova -------------------------
% extract strain and dose from gname
GN = T.(factor1Name);
a = regexpcellout(GN,'_','split');
dose = a(:,2);
dose(cellfun(@isempty,dose)) = {'0mM'};
strain = a(:,1);


T2 = T;
T2.(factor1Name) = [];
DVnames = T2.Properties.VariableNames;
%%
dvn = DVnames';
a = regexprep(dvn,rmfactorName,'');
b = cellfun(@str2num,a);
dmFactorNum = b;

%%
TM = table;
TM.(char(factor1Name)) = GN;
TM.strain = strain;
TM.dose = dose;
IVnames = TM.Properties.VariableNames;

TM = [TM T2];


%% pair wise table
dmTable = table(dmFactorNum,'VariableNames',{rmfactorName});


