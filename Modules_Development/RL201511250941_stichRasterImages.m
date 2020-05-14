% stich raster images
addpath('/Users/connylin/Dropbox/MATLAB/Function_Library_Public');
pHome = '/Users/connylin/Dropbox/RL/MWT_Analysis_ByExp/slo1 60ISI';
[~,~,~,pG] = dircontent(pHome);
pG(regexpcellout(pG,'graffle')) = [];
for gi = 1%:numel(pG)
    p = sprintf('%s/rasterPlot_colorSpeed',pG{gi});
    [~,p] = dircontent(p,'rasterPlot*');
    %%
    for pi = 1%:numel(p)
        load(p{pi})
    end
    
end



