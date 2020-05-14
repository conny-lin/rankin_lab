D = MWTSet.Data_Plate;
D = innerjoin(MWTSet.Info.MWTDB,D);
gu =unique(D.groupname_short);
i = regexpcellout(gu,'N2');
gu = [gu(i);gu(~i)];
su = unique(D.strain);
i = regexpcellout(su,'N2');
su = [su(i);su(~i)];
cu = unique(D.condition_short);
markercolor = nan(size(D,1),3);
cdc = [0 0 0; 1 0 0; 0.5 0.5 0.5; 1 0 1];
for cui = 1:numel(cu)
    i = ismember(D.condition_short,cu(cui));
    markercolor(i,:) = repmat(cdc(cui,:),sum(i),1);
end
close;
fig1 = figure('visible','off'); hold all;
markertype = {'o','^'};
for si = 1:numel(su)
    i = ismember(D.strain,su(si));
    x = D.curve(i);
    y = D.speed(i);
    a = repmat(30,numel(x),1);
    c = markercolor(i,:);
    if si==1
        scatter(x,y,a,c,markertype{si},'filled')
    else
        scatter(x,y,a,c,markertype{si})
    end
    hold on;
end
% label
xlabel('Curvataure'); ylabel('Speed');
% calculate outlier
T = outlierlim(D.curve,D.groupname_short,'multiplier',2);
T.Properties.VariableNames(1) = {'groupname_short'};
T.Properties.VariableNames(2) = {'upperlimx'};
T.Properties.VariableNames(3) = {'lowerlimx'};
T1 = outlierlim(D.speed,D.groupname_short,'multiplier',2);
T1.Properties.VariableNames(1) = {'groupname_short'};
T1.Properties.VariableNames(2) = {'upperlimy'};
T1.Properties.VariableNames(3) = {'lowerlimy'};
A = innerjoin(T,T1);
recx = A.lowerlimx;
recy = A.lowerlimy;
recw = A.upperlimx - A.lowerlimx;
rech = A.upperlimy - A.lowerlimy;
data = [recx recy recw rech];
gname = A.groupname_short;
linestyles = {'-','--'};
for si = 1:numel(su)
    for ci = 1:numel(cu)
        i = ismember(gname,strjoin([su(si) cu(ci)]));
        rectangle('Position',data(i,:),'EdgeColor',cdc(ci,:),...
            'LineStyle',linestyles{si})
    end
end
savefigpdf('speed x curve scatter',pSave);