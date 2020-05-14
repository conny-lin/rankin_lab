msrlist = fieldnames(MWTSet.Data_Group);
D = MWTSet.Data_Group.(msrlist{1});
strains = unique(D.strains);
rowsname = 'predose';
colsname = 'postdose';

rowvaru = unique(D.(rowsname));
colvaru = unique(D.(colsname));

X = repmat((0:0.1:(0.1*(numel(rowvaru)-1)))', 1,numel(colvaru));
yEmpty = nan(numel(rowvaru),numel(colvaru));

% color: black and red
color = [0 0 0; [0.635294139385223 0.0784313753247261 0.184313729405403]];
for msri = 1:numel(msrlist)
    msr = msrlist{msri};
    D = MWTSet.Data_Group.(msr);

for si = 1:numel(strains) % run through each strain
    sname = strains{si};
    % rearrange to output
    Y = yEmpty;
    E = Y;
    for ri = 1:numel(rowvaru)
    for ci = 1:numel(colvaru)
        i = ismember(D.strains, sname) ...
            & ismember(D.(rowsname), rowvaru{ri}) ... % get predose
            & ismember(D.(colsname), colvaru{ci});
        if sum(i)==1; 
            Y(ri,ci) = D.mean(i);
            E(ri,ci) = D.se(i);
        end
        
    end
    end
    
    % create figure
    close;
    f1 = figure('Visible','off');
    axes1 = axes('Parent',f1,'Box','off','FontSize',10);
    hold on;
    ebar1 = errorbar(X,Y,E,'LineStyle','none','Marker','o','MarkerSize',6);
    for ei = 1:size(X,1) 
        c = color(ei,:);
        set(ebar1(ei),'Color',c,'DisplayName',colvaru{ei})
    end
    set(axes1,'XTick',X(:,1),'XTickLabel',rowvaru);
    ylabel(msr);
    xlim([-0.05 X(end)+0.01])
    lg1 = legend('show',colvaru);
    set(lg1,'EdgeColor','none','Location','northeastoutside')
    title(sname);
    % save fig
    savename = sprintf('%s %s',msr,sname);
    printfig(savename,pSave,'closefig',1,'w',3,'h',3)
    
end
end