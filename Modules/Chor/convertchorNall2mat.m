function [pFiles,pMWTval] = convertchorNall2mat(pMWT,inputname,outputname)
%%
% example:
%     outputname = 'Gangnam';
%     inputname = 'gangnam';

%% convert gangnam into .mat
[~,~,pMWTmiss] = getpath2chorfile(pMWT,[outputname,'.mat']);
[pFiles,pMWTval] = getpath2chorfile(pMWTmiss,['*.',inputname,'.*.dat']);
%%

nfiles = numel(pFiles);
% search for files

for mwti = 1:nfiles
    if ismember(mwti,[1:10:nfiles nfiles]) == 1
        fprintf('%d/%d MWT file\n',mwti,nfiles);
    end
    pmwt = pMWTval{mwti};
    pF = pFiles{mwti};
    if ischar(pF); pF = {pF}; end
    nRow = numel(pF);
    
    Data = cell(nRow,1);
    time = nan(nRow,2);
    for wrmi = 1:numel(pF)
        d = dlmread(pF{wrmi});
        Data{wrmi} = d;
        time(wrmi,1:2) = [d(1,1) d(end,1)];
    end
    
    cd(pmwt);
    save([outputname,'.mat'],'Data','time');
end