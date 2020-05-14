
%% transfer files
pDataSource = '/Volumes/Public/Data/MWT_20151111/ExpZip';
pDataBase = '/Volumes/COBOLT/MWT';
% pDataSource = '/Volumes/IRONMAN/MWT_test1';
% pDataBase = '/Volumes/IRONMAN/MWT_test2';
pZip = '/Volumes/IRONMAN/MWT_ZipTemp';
if isdir(pZip) == 0; mkdir(pZip); end
xtart = 3067;

%% process
addpath('/Users/connylin/Dropbox/MATLAB/Function_Library_Public');
addpath('/Users/connylin/Dropbox/RL/Code/Modules/MWTDatabase');
% tsfNewFiles2SameMWTFolder(pDataSource,pDataBase);

%% get data home
addpath('/Users/connylin/Dropbox/RL/Code/Modules/MWTDatabase');
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




