function [pMWT,fMWT] = getpMWT(p,displayoption)
% p = pDataSource;
addpath('/Users/connylin/Dropbox/MATLAB/Function_Library_Public');

if nargin == 1;
    displayoption = 0;
end
if displayoption == 1
    fprintf('getting in from new MWT folder...');
end
[~,~,~,pEf] = dircontent(p);
pEf(regexpcellout(pEf,'__MACOSX')) = [];
% get group folder
[~,~,~,pG] = cellfun(@dircontent,pEf,'UniformOutput',0);
pG(cellfun(@isempty,pG)) = [];
pG = celltakeout(pG,'multirow');
% check if any group folders are actually mwt file/folder
[~,f] = cellfun(@fileparts,pG,'UniformOutput',0);
i = regexpcellout(f,'(\<\d{8}[_]\d{6})|(\<\d{8}[_]\d{6}[.]zip\>)');
pMWTD = pG(i);
pG(i) = [];
[~,~,~,pMWT] = cellfun(@dircontent,pG,'UniformOutput',0);
pMWT(cellfun(@isempty,pMWT)) = [];
pMWT = celltakeout(pMWT,'multirow');
[~,fMWT] = cellfun(@fileparts,pMWT,'UniformOutput',0);
i = regexpcellout(fMWT,'(\<\d{8}[_]\d{6}\>)|(\<\d{8}[_]\d{6}[.]zip\>)');
pMWT(~i) = [];
pMWT = [pMWTD;pMWT];
[~,fMWT] = cellfun(@fileparts,pMWT,'UniformOutput',0);
if displayoption == 1
    fprintf('\n');
    fprintf('-%d mwt files found\n',numel(pMWT));
end
