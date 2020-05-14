function [pFiles,pMWTval] = convertchor2mat(pMWTS,legend,inputname,outputname)
%%
% example:
%     outputname = 'Gangnam';
%     inputname = 'gangnam';

%% convert gangnam into .mat
[~,~,pMWTmiss] = getpath2chorfile(pMWTS,[outputname,'.mat']);
[pFiles,pMWTval] = getpath2chorfile(pMWTmiss,['*',inputname,'.dat']);
nfiles = numel(pFiles);
% search for files
for mwti = 1:nfiles
    if ismember(mwti,[1:10:nfiles nfiles]) == 1
        fprintf('%d/%d MWT file\n',mwti,nfiles);
    end
    pmwt = pMWTval{mwti};
    pF = pFiles{mwti};
    d = dlmread(pF);
    if size(d,2) ~= numel(legend)
        error('data col size does not match legend');
    end
    Data = d;
    time = [d(1,1) d(end,1)];
    cd(pmwt);
    save([outputname,'.mat'],'Data','time','legend');
end