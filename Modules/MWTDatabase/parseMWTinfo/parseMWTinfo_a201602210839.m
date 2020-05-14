function T = parseMWTinfo(pMWTD)
%% accmodate for entry sizes
if iscell(pMWTD) == 0 && numel(pMWTD) == 1
    pMWTD = {pMWTD};
end
if iscell(pMWTD) == 1 && size(pMWTD,1) == 1 && numel(pMWTD) > 1
    pMWTD = pMWTD';
end


%% parse pMWT
if iscell(pMWTD) == 1
    [pG,fMWT] = cellfun(@fileparts,pMWTD,'UniformOutput',0);
    [pE,fG] = cellfun(@fileparts,pG,'UniformOutput',0);
    [~,fE] = cellfun(@fileparts,pE,'UniformOutput',0);
elseif ischar(pMWTD) == 1
    [pG,fMWT] = fileparts(pMWTD);
    fMWT = {fMWT};
    pMWTD = {pMWTD};
    [pE,fG] = fileparts(pG);
    fG = {fG};
    pG = {pG};
    [~,fE] = fileparts(pE);
    fE = {fE};
    pE = {pE};
end


% MWT index
T = table;
T.mwtid = [1:numel(fMWT)]';
T.mwtname = fMWT;
T.mwtpath = pMWTD;
T.expname = fE;
T.exp_date = cellfun(@str2num,regexpcellout(fE,'\<\d{8}','match'));
T.tracker = char(regexpcellout(fE,'(?<=\<\d{8})[A-Z]','match'));
T.expter = char(regexpcellout(fE,'(?<=\<\d{8}[A-Z][_])[A-Z]{2}','match'));
T.groupname = fG;
T.strain = regexpcellout(fG,'\<[A-Z]{1,}\d{1,}','match');
a = regexpcellout(fG,'_','split');
T.rx = regexpcellout(fG,'(?<=\<[A-Z]{1,}\d{1,}_)\w{1,}','match');
T.rx(cellfun(@isempty,T.rx)) = {'NA'};
a = regexpcellout(fE,'_','split');
T.rc = a(:,3);
a = regexpcellout(T.rc,'\d{1,}(?=s)','match');
T.preplate = cellfun(@str2num,a(:,1));
T.ISI = cellfun(@str2num,a(:,2));
T.postrec = cellfun(@str2num,a(:,3));
T.tapN = cellfun(@str2num,regexpcellout(T.rc,'\d{1,}(?=x)','match'));
