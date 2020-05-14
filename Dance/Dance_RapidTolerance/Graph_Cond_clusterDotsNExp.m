close;
postfixname = 'clusterdots Nexp';
D = MWTSet.Data_Exp;
msrlist = {'speed','curve'};
lightblue = [0.301960796117783 0.745098054409027 0.933333337306976];
for msri= 1:numel(msrlist)
    msr = msrlist{msri};
    su = unique(D.(msr).strain);
    for si = 1:numel(su)
        strain = su{si};
        i = ismember(D.(msr).strain,strain);
        D1 = D.(msr)(i,:);
        D1 = sortrows(D1,{'predose','postdose'});
        
        a = strjoinrows(table2cell(D1(:,{'predose','postdose'})),'-');
        gname = regexprep(a,'mM','');
        clusterDotsErrorbar(D1.mean,gname,'markersize',6,...
            'xname','','yname',msr,'scatterdotcolor',lightblue,...
            'visible','off');
        title(strain);
        savename = sprintf('%s %s %s',msr,postfixname,strain);
        printfig(savename,pSave,'w',4,'h',4,'ext','pdf','closefig',0)
    end
end
  