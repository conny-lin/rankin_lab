function figure1 = create_colorlegend(minum,maxnum,colormapname)

%% defaults


%% create color scheme
if nargin<=2
x = [-1:0.01:0 0.01:0.01:1];
color_scheme = repmat(x,50,1);

n = numel(x);
xticks = [1 21:20:n];

n2 = (n-1)/2;
xlblist = [minum:-minum/n2:0 maxnum/n2:maxnum/n2:maxnum];

xlables = xlblist(xticks);
xlables = regexprep(cellstr(num2str(xlables'))',' ','');


figure1 = figure;
axes1 = axes('Parent',figure1,...
    'XTickLabel',xlables,...
    'XTick',xticks,...
    'TickDir','out',...
    'Layer','top', 'fontSize',16);
hold on
imagesc(color_scheme)
ylim([1 50])
xlim([1 201])

end

if nargin >=3
    
    if isstr(colormapname)
       cmap = colormap(colormapname);
    else
        cmap = colormapname;
    end
     
    x = linspace(minum,maxnum,size(cmap,1));
    
    color_scheme = repmat(x,50,1);

    n = numel(x);
    xticks = [1 21:20:n];

    n2 = (n-1)/2;
    xlblist = [minum:-minum/n2:0 maxnum/n2:maxnum/n2:maxnum];

    xlables = xlblist(xticks);
    xlables = regexprep(cellstr(num2str(xlables'))',' ','');


    figure1 = figure;
    axes1 = axes('Parent',figure1,...
        'XTickLabel',xlables,...
        'XTick',xticks,...
        'TickDir','out',...
        'Layer','top', 'fontSize',16);
    hold on
            colormap(colormapname);

    imagesc(color_scheme)
    ylim([1 50])
    xlim([1 201])
    
end