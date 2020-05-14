function [DataS,tptable] = DataPrep_rANOVA(Data,factorlist)

%%
if nargin==1
factorlist = 'tap~mwtid@RevFreq+RevSpeed+RevDur';
end
% vararginProcessor

%% translate factorlist
a = regexp(factorlist,'~','split');
rmfactor = a{1};
a = regexp(a{:,2},'@','split');
gfactor = a{1};
msrlist = regexp(a{:,2},'+','split');


%%
tpu = unique(Data.(rmfactor));
measNames = num2cellstr(tpu);
tpnames = strjoinrows([cellfunexpr(measNames,rmfactor),measNames],'');


mwtidu = unique(Data.(gfactor));
A = nan(numel(mwtidu)*numel(msrlist),numel(tpu));
Leg = cell(numel(mwtidu)*numel(msrlist),2);
row1 = 1;
for mwti = 1:numel(mwtidu)
    for msri = 1:numel(msrlist)
        d = table2array(Data(Data.mwtid==mwtidu(mwti),msrlist(msri)));
        col = table2array(Data(Data.mwtid==mwtidu(mwti),{'tap'}));
        A(row1,col) = d;
        Leg{row1,1} = msrlist{msri};
        Leg{row1,2} = mwtidu(mwti);
        row1 = row1+1;
    end
end

%%
A = array2table(A,'VariableNames',tpnames);
Leg = cell2table(Leg,'VariableNames',{'msr',gfactor});

DataS = [Leg A];
tptable = table(tpu,'VariableNames',{'tap'});

