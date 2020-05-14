% get Dance_Gee graph out
pHome = '/Users/connylin/Dropbox/RL/MWT_Analysis_ByExp/60sISI by Strain/Strains';
pNew = [fileparts(pHome),'/StrainGleeGraphs'];
if isdir(pNew)==0; mkdir(pNew); end

strainlist = dircontent(pHome);

for x = 1:numel(strainlist)
    strain = strainlist{x};
    p = sprintf('%s/%s/Dance_Glee/Graph HabCurve/%s',pHome,strain,strain);
    [f2,p2] = dircontent(p);
    for y = 1:numel(p2)
        ps = p2{y};
        pd = sprintf('%s/%s',pNew,f2{y});
        copyfile(ps,pd,'f');
    end
end

