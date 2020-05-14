function ephys_graphline(DataG,pSave,tap1,tap2,strain)

%% graph space out graphs
% decide group sequence
gu = cell(size(DataG,2),1);
for gi =1:size(DataG,2)
    gu{gi} = DataG(gi).name;
end
guseq = [find(regexpcellout(gu,'N2'));find(~regexpcellout(gu,'N2'))]';


% graph setting
linespec = {'.k','.r','-k','-r'};
linespec = linespec(guseq);
% graph mean curve 
close
figure('Visible','off'); hold on;
for gi = guseq
    gn = regexprep(DataG(gi).name,'_',' ');
    d = DataG(gi).speedbm;
    x = DataG(gi).time(1,:);
    y = nanmean(d);
    y(x==0) = nan; % make t=0 nan
%     sd = nanstd(d);
%     n = sum(~isnan(d));
%     se = sd./sqrt(n-1);
    e1 = plot(x,y,linespec{gi},'DisplayName', gn,'LineWidth',1,'MarkerSize',1);        
end

xlim([-.5 2]);
ylim([-.3 .4]);
ylabel('velocity (body length/s)')
title(strain);
% legend1 = legend('show');
% set(legend1,'Location','southeast','EdgeColor',[1 1 1]);
savename = sprintf('%s ephys t%d-%d pubgrade',strain,tap1,tap2);

printfig(savename,pSave,'w',2,'h',2,'closefig',0);