function fig = graphLine(msrlist)

%% graph ===============================
for msri = 1:numel(msrlist)

    msr = msrlist{msri};
    
    % get graph data ----------------
    D = G.(msr);
    y = D.Y;
    e = D.E;
    x = D.X;
    gn = D.groupname;
    D1(:,:,1) = x;
    D1(:,:,2) = y;
    D1(:,:,3) = e;
    % sort by N2
    gns = sortN2first(gn,gn);
    [i,j] = ismember(gn,gns);
    D1 = D1(:,j,:);
    X = D1(:,:,1);
    Y = D1(:,:,2);
    E = D1(:,:,3);
    % -------------------------------

    % setting --------------------------
    w = 4.5;
    h = 3.5;
    gp = graphsetpack('cathyline');
    gnss = regexprep(gns','_',' ');
    gnss = regexprep(gnss,strain,genotype);
    gnss = regexprep(gnss,'N2','Wildtype');
    gp.DisplayName = gnss;
    gpn = fieldnames(gp);
    % -----------------------------------

    % plot ------------------------------------------
    fig = figure('Visible','off','PaperSize',[w h],'Unit','Inches',...
        'Position',[0 0 w h]); 
    ax1 = axes('Parent',fig,'XTick',[0:10:30]); 
    hold(ax1,'on');
    e1 = errorbar(X,Y,E,'Marker','o','MarkerSize',4.5);
    % -----------------------------------------------

    % axis +++++++++
    axs = get(ax1);
    ax1.TickLength = [0 0];
    ax1.FontSize = 10;
    % -------------

    % get settings -----------------------------------
    for gi = 1:numel(gn)
        for gpi = 1:numel(gpn)
            e1(gi).(gpn{gpi}) = gp.(gpn{gpi}){gi};
        end
    end
    % -----------------------------------------

    % title -------------------------
    title(genotype)
    % -------------------------------

    % x axis ++++++++
    xlabel('Tap')
    xlim([0 30.5]);
    % --------------

    % y axis +++++++
    yname = sprintf('P (%s)',rTargetName);
    ylim([0 1]);
    ylabel(yname);
    % ----------------


end
% =========================================================================
