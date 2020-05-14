% get Dance_Gee graph out
pDes = '/Users/connylin/Dropbox/RL/MWT_Analysis_ByExp/STH/10sISI by Strain/Strains';
[strainlist, pStrain] = dircontent(pDes);

pHome = '/Users/connylin/Dropbox/RL/MWT_Analysis_ByExp/STH/by pathway';
[~,~,~,pFF] = dircontent(pHome);
ignoreList = {'DA465'};
for y = 1:numel(pFF)
    pFig = pFF{y};
    [~,fn] = fileparts(pFig);
    fprintf('\nProcessing [%s]\n',fn);
    [~,~,~,p2] = dircontent(pFig);
%     [~,~,~,p] = cellfun(@dircontent,p,'UniformOutput',0);
%     p = celltakeout(p);
    i = regexpcellout(p2,'(graffle)|(N2)|(N2_400mM)|(^Dance)');
    p2(i) = [];
    if numel(a) < 2
       p2
       error('might not have 400mM group'); 
    end
    
    for x =1:numel(p2)
        
       p = p2{x};
       [~,gn] = fileparts(p);
       if sum(ismember(gn,ignoreList)) == 0
           fprintf('-[%s]\n',gn);
           cond = ~isempty(regexp(gn,'400mM'));
           if cond == 1
               strain = regexprep(gn,'_400mM','');
               cond = '400mM';
           else
               strain = gn;
               cond = '0mM';
           end
           i = ismember(strainlist,strain);
           if sum(i) ~= 1
              warning('no home strain folder found'); 
           else

               pdH = sprintf('%s/rasterPlot_colorSpeed/%s %s',pStrain{i},strain,cond);   
               pF = [p,'/rasterPlot_colorSpeed'];
               if isdir(pF) == 0;
                   x
                   error('no raster plot');
               end
               [~,ps] = dircontent(pF);
               n = numel(ps);
               pd = regexprep(ps,pF,pdH);
               p = char(unique(cellfun(@fileparts,pd,'UniformOutput',0)));
               if isdir(p) == 0; mkdir(p); end

               cellfun(@copyfile,ps,pd);
           end
       end
    end
end

