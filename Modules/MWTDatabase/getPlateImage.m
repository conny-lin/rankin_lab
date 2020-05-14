function getPlateImage(pMWT,varargin)


%% DEFAULTS AND PROCESS VARARGIN
pSave = '/Users/connylin/Dropbox/RL/Dance Output'; % default save output path
% process varargin
vararginProcessor;

%% get photo
for mi = 1:numel(pMWT)
    pmwt = pMWT{mi};
    [f,p] = dircontent(pmwt,'*.png');
    try
        copyimage(p,pSave)
    catch err 
        p(regexpcellout(f,'^[.]'))= [];  % ignore temp files
        if numel(p) > 1
            disp(char(p))
            [fset,pset] = dircontent(pmwt,'*.set');
            if numel(fset) > 1
                rethrow(err)
            else
                p = regexprep(pset,'[.]set\>','.png');
                copyimage(p,pSave)
            end
        end
    end
end

end

function copyimage(p,pSave)
    ps = char(p);
    a = parseMWTinfo({fileparts(ps)});
    pd = sprintf('%s/%s E%d %s.png',pSave,char(a.groupname), a.exp_date, char(a.mwtname));
    copyfile(ps,pd); 
end