close;
D = MWTSet.Data_Plate;
D = innerjoin(MWTSet.Info.MWTDB,D);
D.exp_datestr = cellstr(num2str(D.exp_date));
msrlist = {'speed','curve'};
for msri= 1:numel(msrlist)
    msr = msrlist{msri};
    D1 = D.(msr);
    gname = strjoinrows(D(:,{'exp_datestr','groupname_short'}));
    clusterDotsErrorbar(D1,gname,'markersize',5,'xname','','yname',msr);
    xticklabel_rotate([],90);
    savefigpdf(sprintf('%s bar scatter by plate sep exp',msr),pSave);    
end