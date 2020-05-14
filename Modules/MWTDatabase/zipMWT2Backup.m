function zipMWT2Backup(pDataBase,pZip,varargin)
%% zipMWT2Backup(pDataSource,pDataBase) 
% make MWT zip file contianing only blob, png, set and summary files to
% back up folder
% 
% Input:
%     pDataBase = '/Volumes/COBOLT/MWT';
%     pZip = '/Volumes/IRONMAN/MWT_ZipTemp';

%% function paths
addpath('/Users/connylin/Dropbox/MATLAB/Function_Library_Public');
addpath('/Users/connylin/Dropbox/RL/Code/Modules/MWTDatabase');
% tsfNewFiles2SameMWTFolder(pDataSource,pDataBase);

%% DEFAULTS, VARARGIN, AND OTHER SETTINGS
% defaults
nInput = 2;
xtart = 1;

% varargin processer
if nargin > nInput
    A = varargin;
    if isinteger(numel(A)/2) == 1
        error('variable entries incorrect')
    end
    callName = A(1:2:numel(A));
    % check if call names are strings
    for x = 1:numel(callName)
       if ischar(callName{x}) == 0
           error('varargin incorrect');
       end
    end
    % get values
    setValue = A(2:2:numel(A));
    for x = 1:numel(callName)
        % determine input type
        if eval(sprintf('ischar(%s)',callName{x})) == 1
            eval(sprintf('%s = ''%s'';',callName{x},setValue{x}));
        elseif eval(sprintf('isnumeric(%s)',callName{x}))
            eval(sprintf('%s = %d;',callName{x},cell2mat(setValue(x))));
        elseif eval(sprintf('iscell(%s)',callName{x}))
            a = setValue{x};
            eval(sprintf('%s = a;',callName{x}));
        else
            error('need coding');
        end
    end
end

% process inputs
if isdir(pZip) == 0; mkdir(pZip); end



%% get data home
[pMWT,fMWT] = getpMWT(pDataBase);
% check if duplicated files in database folder
a = tabulate(fMWT);
i = cell2mat(a(:,2))<=1;
a(i,:) = [];
if isempty(a) == 0
    disp(a);
    error('duplicate MWT files');
end

%% zip MWT files
reportJump = [1 100:100:numel(pMWT)];
for x = xtart:numel(pMWT)
    if sum(reportJump == x) > 0
        fprintf('Processing %d/%d MWT files\n',x,numel(pMWT));
    end
    pmwt = pMWT{x};
    [~,fmwt] = fileparts(pmwt);
    % check if zip already exist
    pd = regexprep(pmwt,pDataBase,pZip);
    pdH = fileparts(pd);
    if isdir(pdH) == 0; 
        mkdir(pdH); 
        a = {};
    else
        a = dircontent(pdH,[fmwt,'.zip']);
    end
    if isempty(a) == 1
        fprintf('-zipping [%s]  ',fmwt);
        [mwtcontent,pmwtC] = dircontent(pmwt);
        i = regexpcellout(mwtcontent,'([.]blob)|([.]blobs)|([.]png)|([.]set)|([.]summary)');
        if sum(i) > 0
            fprintf('%d files  ',sum(i));
            ps = pmwtC(i);
            if sum(i) == 1
                psf = fileparts(char(ps));
                zip(pd,psf);
            else
                zip(pd,ps);
            end
            fprintf('done\n');
        else
            fprintf('0 files can not zip\n');
        end
    else
        fprintf('-[%s] already zipped\n',fmwt)
    end
    
end




