function zipMWTclean(pMWTnew,pStd,pCheck)
%%
% pMWTnew - pmwt 
% p standard - pdata home
% pcheck - zip file designation

%%
for y = 1:numel(pMWTnew)
    fprintf('zipping; %d/%d\n',y,numel(pMWTnew));

    %% get path
    pmwt = pMWTnew{y};
    
    %% check integrity of MWT files
    [~,pfiles] = chekMWTintegrity(pmwt);
        
    %% zip to designation
    pd = regexprep(pmwt,pStd,pCheck);
    if isdir(fileparts(pd))==0; mkdir(fileparts(pd)); end
    zip([pd,'.zip'],pfiles)
end