function txt = reportTAR(A)

%%  variables =============================================================
pvsig = 0.05;
pvlim = 0.001;
rspName = 'acceleration response probability';
controlname = 'wild-type';
% get genotype
strain = char(A.strain);
strainNames = DanceM_load_strainInfo(strain);
genotype = strainNames.genotype{ismember(strainNames.strain, strain)};
% =========================================================================



%% anova sig % =============================================================
txt = 'A repeated measures ANOVA evaluating the effects of strain, and ethanol on acceleration response probability revealed';
D = A.anova;
pv = D.pvalue;
i = pv < pvsig;
if sum(i) >= 1
    F = [D.factor(i) D.F(i) D.ptxt(i)];
    if size(F,1) == 1
       Fsigtxt = strjoin(F',', ');
    else
        a = F(1:end-1,:);
        b = strjoinrows(a,', ');
        c = strjoin(b',', ');
        d = sprintf(', and %s',strjoin(F(end,:),', '));
        Fsigtxt = [c d];
    end
    s = sprintf('significant main effects of %s',Fsigtxt);
    s = regexprep(s,'dose','ethanol');
    txt = sprintf('%s %s',txt,s);
    siganova = true;

else
    siganova = false;
end

% anova not sig
pv = D.pvalue;
i = pv > pvsig;
if sum(i) > 1
    F = [D.factor(i)];
    if size(F,1) == 1
       Fsigtxt = strjoin(F',', ');
    else
        a = F(1:end-1,:);
        b = strjoinrows(a,', ');
        c = strjoin(b',', ');
        d = sprintf(' or %s',strjoin(F(end,:),', '));
        Fsigtxt = [c d];
    end
    if siganova
        s = sprintf('but not %s',Fsigtxt);
    else
        s = sprintf('no significant effect of %s',Fsigtxt);
    end
    txt = sprintf('%s, %s',txt,s);
end
txt = sprintf('%s.',txt);
% =========================================================================



%% comparison between mutant 0mM and widltype 0mM =========================
txt = sprintf('%s In the absence of ethanol, the %s mutant',txt,genotype);
D = A.posthoc_groups;
compname = sprintf('N2*%s', strain);
pv = D.pvalue(compname);
pvt = char(D.ptxt(compname));
% per tap comparison
a = regexpcellout(compname,'*','split');
posthocstr = genPosthocTapbyGStr(A,a{1},a{2},pvlim,pvsig);

if pv < pvsig
    s = sprintf('had different %s compared to the %s 0mM group',rspName,controlname);
    s = sprintf('%s (curve, %s, pairwise comparison, %s)',s,pvt,posthocstr);
    wt0mt0sig = true;
else
    s = sprintf('had similar %s compared to the %s 0mM group',rspName,controlname);
    wt0mt0sig = false;

end
txt = sprintf('%s %s.',txt,s);
% =========================================================================



%% comparison in wildtype =================================================
txt = sprintf('%s In the presence of ethanol, %s',txt,controlname);
D = A.posthoc_groups;
compname = {'N2*N2_400mM'};
pv = D.pvalue(compname);
pvt = char(D.ptxt(compname));
% per tap comparison
a = regexpcellout(compname,'*','split');
posthocstr = genPosthocTapbyGStr(A,a{1},a{2},pvlim,pvsig);

if pv < pvsig
    s = sprintf('showed the expected elevation of %s',rspName);
    s1 = sprintf('(curve, %s, pairwise comparison, %s)',pvt,posthocstr);
    s = sprintf('%s %s',s,s1);
    wtsig = true;
else
    s = sprintf('failed to show the expected elevation of %s',rspName);
    wtsig = false;
end

txt = sprintf('%s %s.',txt,s);
% =========================================================================



%% comparison in mutants =================================================
txt = sprintf('%s Comparison between the %s 0mM and 400mM groups',txt,genotype);
D = A.posthoc_groups;
compname = sprintf('%s*%s_400mM',strain, strain);
pv = D.pvalue(compname);
pvt = char(D.ptxt(compname));
% per tap comparison
a = regexpcellout(compname,'*','split');
posthocstr = genPosthocTapbyGStr(A,a{1},a{2},pvlim,pvsig);
clear s;
if pv < pvsig
%     if wtsig
%         s = sprintf('and the');
%     else
%         s = sprintf('but did for the');
%     end
    s = sprintf('showed that ethanol altered %s (curve, %s, pairwise comparison, %s)',rspName, pvt,posthocstr);
    mutsig = true;
else
%     if wtsig
%         s = sprintf('but ethanol had no effect on the %s mutant', genotype);
%     else
%         s = sprintf('nor did the %s mutant', genotype);
%     end
    s = sprintf('showed that ethanol did not alter %s', rspName);

    mutsig = false;
end
txt = sprintf('%s %s.',txt,s);
% =========================================================================



%% comparison between mutant 400mM and widltype 400mM =====================
if mutsig 
    D = A.posthoc_groups;
    compname = sprintf('N2_400mM*%s_400mM', strain);
    pv = D.pvalue(compname);
    pvt = char(D.ptxt(compname));
    % per tap comparison
    a = regexpcellout(compname,'*','split');
    posthocstr = genPosthocTapbyGStr(A,a{1},a{2},pvlim,pvsig);
    % constrct text
    s = sprintf('Comparison between the %s 400mM and the %s 400mM group showed that',controlname,genotype);
    if pv < pvsig
        s = sprintf('%s the %s 400mM group had a different %s',s, genotype,rspName);
        s = sprintf('%s (curve, %s, pairwise comparison, %s)',s,pvt,posthocstr);
    else
        s = sprintf('%s %s had similar %s',s, genotype, rspName);
    end
    txt = sprintf('%s %s.',txt,s);
end
% =========================================================================



%% comparison between mutant 400mM and widltype 0mM =======================
if ~mutsig 
    D = A.posthoc_groups;
    compname = sprintf('N2*%s_400mM', strain);
    if sum(ismember(D.g1,'N2'))== 0
        compname = sprintf('%s_400mM*N2', strain);
    end
    pv = D.pvalue(compname);
    pvt = char(D.ptxt(compname));
    % per tap comparison
    a = regexpcellout(compname,'*','split');
    posthocstr = genPosthocTapbyGStr(A,a{1},a{2},pvlim,pvsig);
    % constrct text
    s = sprintf('Comparing to the %s 0mM group, %s on ethanol',controlname,genotype);
    if pv < pvsig
        s = sprintf('%s had different %s',s,rspName);
        s = sprintf('%s (curve, %s, pairwise comparison, %s)',s,pvt,posthocstr);
    else
        s = sprintf('%s %s had similar %s',s, genotype, rspName);
    end
    txt = sprintf('%s %s.',txt,s);
end
% =========================================================================



%% comparison between mutant 0mM and widltype 400mM =======================
if ~mutsig && wtsig && wt0mt0sig
    D = A.posthoc_groups;
    compname = sprintf('N2_400mM*%s', strain);
    if sum(ismember(D.g1,'N2_400mM'))== 0
        compname = sprintf('%s*N2_400mM', strain);
    end
    pv = D.pvalue(compname);
    pvt = char(D.ptxt(compname));
    % per tap comparison
    a = regexpcellout(compname,'*','split');
    posthocstr = genPosthocTapbyGStr(A,a{1},a{2},pvlim,pvsig);
    % constrct text
    s = sprintf('In fact, comparing to the %s 0mM group, the %s 0mM group',controlname,genotype);
    if pv < pvsig
        s = sprintf('%s had different %s',s,rspName);
        s = sprintf('%s (curve, %s, pairwise comparison, %s)',s,pvt,posthocstr);
    else
        s = sprintf('%s %s had similar %s',s, genotype, rspName);
    end
    txt = sprintf('%s %s.',txt,s);
end
% =========================================================================



%% conclusion =============================================================
if ~mutsig && ~wtsig
   s = sprintf('Because the %s mutant did not show an effect of ethanol, the %s mutation eliminated ethanol''s effect on %s',genotype, genotype,rspName);
    txt = sprintf('%s %s.',txt,s);
elseif mutsig && wtsig
    s = sprintf('Because the %s mutant showed an effect of ethanol, the %s mutation is not important for ethanol''s effect on %s',genotype, genotype,rspName);
    txt = sprintf('%s %s.',txt,s);
end
% =========================================================================





























