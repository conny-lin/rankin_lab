function [pSPS,pMWT_success,pMWT_failed] = getpath2chorfile(pMWT,chorfilename,varargin)

%% get paths to chor file with extension name chorfilename
% chorfilename = '.drunkposture2.dat';
% OUTPUT:
%     pSP only outputs valid paths, 

%% examine varargin and defaults
iterationgap = 50;
reporting = true;
vararginProcessor;

%% calculate variables
nfiles = numel(pMWT);
pSP = cell(size(pMWT));
iterationgaplist = [1:iterationgap:nfiles nfiles];

%% search for files

for mwti = 1:nfiles
    % report every 10 iterations
    if ismember(mwti,iterationgaplist) == 1 && reporting
        fprintf('%d/%d MWT file\n',mwti,nfiles);
    end
    pmwt = pMWT{mwti}; % get current folder path
    %% get extenion name
    if strfind(chorfilename,'*')
        fname = chorfilename;
    else
        fname = ['*',chorfilename]; 
    end
    % search for 
    [fn,p] = dircontent(pmwt,fname); 
    i = regexpcellout(fn,'\<[.]');
    if sum(i) > 0
        delete(p{i});
        p(i) = [];
    end
    if isempty(p) == 0
        if numel(p) == 1
            pSP(mwti) = p; % store
        else
            pSP{mwti} = p;
        end
    end
end

%% see which one has files
i = ~cellfun(@isempty,pSP);
pSPS = pSP(i);
pMWT_success = pMWT(i);
pMWT_failed = pMWT(~i);
if reporting
    fprintf('%d/%d mwt files contains %s\n',numel(pMWT_success),nfiles,chorfilename);
end


end


