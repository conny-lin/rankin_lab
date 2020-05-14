% function MWT_fileList = makeMWTchoroutputlist(databasefile,pSave)

%% remake MWTDb
addpath('/Users/connylin/Dropbox/RL/Code/Modules/MWTDatabase');
% pData = '/Volumes/COBOLT/MWT';
% pSave = '/Users/connylin/Dropbox/RL/MWTDB';
% databasefile = '/Users/connylin/Dropbox/rl/MWTDB/MWTDB.mat';


%% load database
load(databasefile)
pMWT = MWTDB.text.mwtpath;


%% survey MWT files
Prefix = table;
Prefix.mwtpath = pMWT;
Prefix.name = cell(size(pMWT));
Prefix.issue = cell(size(pMWT));
FileTemp = cell(numel(pMWT),2);
itrgap = 50;
nfiles = numel(pMWT);
%%
for mwti = 2110:nfiles
    if ismember(mwti,1:itrgap:nfiles) == 1
        fprintf('%d/%d MWT file\n',mwti,nfiles);
    end
    pmwt = pMWT{mwti};
    fn = dircontent(pmwt);
    % get file extensions
    dotparts = regexpcellout(fn,'[.]','split');
    
    % -delete temp files
    i= cellfun(@isempty,dotparts(:,1));
    if sum(i) > 0
        fnT = fn(i);
        % --make sure it's dot at the begining
        i = regexpcellout(fnT,'\<[.]');
        fnT(~i) = [];
        % --delete 
        disp(char(fnT))
%         fprintf('mwt:%d going to delete, press enter to continue...',mwti);
%         pause;
%         fprintf('\n');
        cd(pmwt); cellfun(@delete,fnT);
        fn = fn(i);
    end
    
    % -get suffix
    dotparts = regexpcellout(fn,'[.]','split');
    suffix = cell(size(dotparts,1),1);
    for fi = 1:numel(fn)
        f = dotparts(fi,2:end);
        f(cellfun(@isempty,f)) = [];
        % eliminate number
        f(regexpcellout(f','\d{5}')) = [];
        if numel(f) > 1
            suffix{fi} = strjoin(f,'.');
        elseif isempty(f) == 1
            suffix{fi} = 'NA';
        else
            suffix(fi) = f;
        end
    end
    
    % summarize suffix
    a = tabulate(f);
    sf = a(:,1);
    n = cell2mat(a(:,2));
    nrow = numel(sf);
    FileTemp{mwti,1} = repmat({pmwt},nrow,1);
    FileTemp{mwti,2} = sf;
    FileTemp{mwti,3} = n;

    
    % -store prefix
    a = fn;
    a(regexpcellout(a,'[.]trig\>')) = [];
    a = regexprep(a,'(_\d{5}k[.]blobs\>)|(_\d{5}[.]blob\>)','');
    a = regexpcellout(a,'[.]','split');
    a = a(:,1);
    prefixU = unique(a);
    if numel(prefixU) == 1
        Prefix.name(mwti) = prefixU;
    elseif numel(prefixU) > 1
        Prefix.name(mwti) = prefixU(1);
        Prefix.issue(mwti) = {'multiple'};
    elseif numel(prefixU) == 0
        error('fix no prefix');
    end
end



%% resconsturct

n = cellfun(@numel,FileTemp(:,1));
nRow = unique(sum(n));
rowind = cumsum(n);
rowind = [1;rowind];
mwtpath = cell(nRow,1);
filename = mwtpath;
filecount = nan(nRow,1);
for mwti = 1:size(FileTemp,1)
    p = FileTemp{mwti,1};
    n = numel(p);
    r1 = rowind(mwti);
    r2 = rowind(mwti+1);
    mwtpath(r1:r2) = p;
    filename(r1:r2) = FileTemp{mwti,2};
    filecount(r1:r2) = FileTemp{mwti,3};
end


%% make table
T = table;
T.mwtpath = mwtpath;
T.extname = filename;
T.count = filecount;
MWT_fileList = T;
cd(pSave);
save('MWT_fileList.mat','MWT_fileList');
writetable(MWT_fileList,'MWT_fileList.txt','Delimiter','tab');














