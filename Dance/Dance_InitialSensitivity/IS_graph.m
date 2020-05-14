function IS_graph(DataMeta,MWTDB,pSave,startTime,endTime)
gns = MWTDB.groupname(DataMeta.mwtid);
gnseq = unique(gns);
i = regexpcellout(gnseq,'N2');
gnseq = [gnseq(i);gnseq(~i)];
MWTSet.Data_GroupByWorm = struct;
MWTSet.Data_GroupByWorm.curve = grpstatsTable(DataMeta.curve, gns,'gnameu',gnseq,'gnameutitle','groupname');
MWTSet.Data_GroupByWorm.speedbm = grpstatsTable(DataMeta.speedbm, gns,'gnameu',gnseq,'gnameutitle','groupname');
MWTSet.Data_GroupByWorm.Speed = grpstatsTable(DataMeta.speed, gns,'gnameu',gnseq,'gnameutitle','groupname');

%% graph
GraphData = MWTSet.Data_GroupByWorm;
msrlist = fieldnames(MWTSet.Data_GroupByWorm);
for msri = 1:numel(msrlist)
    close;
    gn = GraphData.(msrlist{msri}).groupname;
    gn = regexprep(regexprep(gn,'_',' '),'mM','');
    y = GraphData.(msrlist{msri}).mean;
    e = GraphData.(msrlist{msri}).se;
    msr = msrlist{msri};
    savename = sprintf('%s t%d-%d',msr,startTime,endTime);

    % plot
    fig1 = figure('Visible','off');
    ax1 = axes('Parent',fig1,'Box','off','XTick',1:numel(gn),...
        'XTickLabel',gn');
    hold(ax1,'on');

    bar(y,'EdgeColor','none','FaceColor',[0 0.447058826684952 0.74117648601532])
    hold on;
    errorbar([1:numel(gn)]',y,e,'Color',[0 0 0],'LineStyle','none','LineWidth',1.2)

    ylabel(msrlist{msri})
    xlim([.5 numel(gn)+.5])
    printfig(savename,pSave,'w',3.5,'h',2,'closefig',1)
end
