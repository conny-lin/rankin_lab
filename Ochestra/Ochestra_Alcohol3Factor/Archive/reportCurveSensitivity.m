function txt = reportCurveSensitivity(curveData)

%%
pvsig = 0.05;
pvlim = 0.001;
%%

strain = curveData.strain;
genotype = curveData.genotype;
% ------------------------


%% opening sentence +++++++++
txt = 'For ethanol sensitivity measured by';
txt = sprintf('%s the percentage of body curve depressed by ethanol,',txt);
% ---------------------------

% wildtype percentage effect +++++++++++++++++++++++++++++
D = percentChange(curveData);
v1 = D.pct('N2');
txt = sprintf('%s wildtype showed a %.0f%%',txt,abs(v1));
if v1<0
    txt = sprintf('%s decrease',txt);
else
    txt = sprintf('%s increase',txt);
end
% ----------------------------------------------------------

% mutant percentage effect +++++++++++++++++++++++++++++
D = percentChange(curveData);
v2 = D.pct(strain);
if v2<0
    effectword = 'decrease';
else
    effectword = 'increase';
end
txt = sprintf('%s, and the %s strain showed a %.0f%% %s',txt, genotype, abs(v2), effectword);
% -------------------------------------------------------

% differences +++++++++++
d = round(abs(v2-v1));
if d==0
    txt = sprintf('%s, indicating % mutant had the same sensitivity to ethanol''s effect on curve',txt,genotype);
else
    if v2<v1 
        dirtext = 'reduction';
    else
        dirtext = 'increase';
    end
    txt = sprintf('%s, representing a %.0f%% %s in ethanol''s effect on body curve compared to wildtype',txt,d,dirtext);
end
% ---------------------

% anova +++++++++++
A = curveData.manova;
i = A.pvalue < pvsig;
a = strjoin(A.txt',', ');
s = sprintf('(ANOVA, %s',a);
txt = sprintf('%s %s',txt, s);
% ------------------

% pairwise +++++++++++
D = percentChange(curveData);
gn = D.Properties.RowNames;
gn =regexprep(gn,'N2','wildtype');
gn = regexprep(gn,strain,genotype);
b = repmat({'0mM vs 400mM'},numel(gn),1);
c = strjoinrows([gn b], ' ');
d = D.p;
s = strjoin(strjoinrows([c d],', ')',', ');
s = sprintf('pairwise comparison, %s)',s);
txt = sprintf('%s %s.',txt, s);
% ----------------------

%% difference between 0mM and 400mM groups ++++++++++++++
D = curveData.posthoc;
i = ismember(D.g1,'N2 0mM') & ismember(D.g2,[strain,' 0mM']);
if sum(i)~=1
    i = ismember(D.g2,'N2 0mM') & ismember(D.g1,[strain,' 0mM']);
end
pv0 = D.pvalue(i);
pvs0 = D.p{i};
if pv0 < pvsig
    wt0mt0sig = true;
else
    wt0mt0sig = false;
end

i = ismember(D.g1,'N2 400mM') & ismember(D.g2,[strain,' 400mM']);
if sum(i)~=1
    i = ismember(D.g2,'N2 400mM') & ismember(D.g1,[strain,' 400mM']);
end
pv4 = D.pvalue(i);
pvs4 = D.p{i};
if pv4 < pvsig
    wt4mt4sig = true;
else
    wt4mt4sig = false;
end
% compose strings
s1 = sprintf('The results indicated that the %s mutation',genotype);
if wt4mt4sig
    s2 = sprintf('altered ethanol''s effect on body curve');
else
    s2 = sprintf('did not altered ethanol''s effect on body curve');

end

s3 = sprintf('(wildtype 0mM vs %s 0mM, %s',genotype,pvs0);
s4 = sprintf('wildtype 400mM vs %s 400mM, %s)',genotype,pvs4);
txts = sprintf('%s %s %s %s',s1,s2,s3,s4);
txt = sprintf('%s %s.',txt,txts); 




%% -----------------------------------------------------
% =========================================================================




























