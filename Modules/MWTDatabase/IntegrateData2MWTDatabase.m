%% MWTDatabase_v3 
%**************************************************************************
% modified from MWTDatabase_v2
% data must be prepared according to this instruction:
% https://www.evernote.com/l/ADe-7p_cJCJPmae-azp-MFU988pObC454iQ
% file format must be like this: 20130409A_CL_100s30x10s10s_slo1
% 
% pNewData = path of folder that contains new data.
% pDataBase = path of folder that contain the MWT database.
%**************************************************************************
function MWTDB = IntegrateData2MWTDatabase(pNewData,pDataBase)
%**************************************************************************
% INITIALIZING 
%--------------------------------------------------------------------------
% clc; clear; close all;
% addpath('/Users/connylin/Dropbox/Code/Matlab/Library/General');
% pM = setup_std(mfilename('fullpath'),'RL','genSave',true);
% addpath(pM);
%**************************************************************************
% SET UP
%**************************************************************************
if nargin==0
    pNewData = '/Volumes/COBOLT/MWT_New';
    pDataBase = '/Volumes/COBOLT/MWT';
end
%**************************************************************************
% CLEAN UP
%--------------------------------------------------------------------------
% Clean database temp files 
%--------------------------------------------------------------------------
% [~,p] = dircontent(pDataBase);
% takeout_tempfiles(p); % take out temp files
%**************************************************************************
% PROCESS NEW FILES 
%**************************************************************************
% experiment names 
% --------------------------------------------------------------------------
[~,pExp] = dircontent(pNewData); % survey for files
% return if no data found
if isempty(pExp)
    fprintf('no new data found');
   return 
end
[pExp,en] = takeout_tempfiles(pExp); % take out temp files
v = regexpcellout(en,'^\d{8}[A-Z][_][A-Z]{2}[_]\d{1,}s\d{1,}x\d{1,}s\d{1,}s'); % name validator
if sum(v) == numel(en) % if names are valid
    fprintf('Validated: experiment names\n');
else
    fprintf('Error: experiment names below are incorrect: \n');
    disp(char(en(~v)))
    error('Please correct above exp names');
end
%--------------------------------------------------------------------------
% check if experiment already exist in database
%--------------------------------------------------------------------------
[fnD,~] = dircontent(pDataBase);
i = ismember(en,fnD);
if any(i)
    warning('some experiments already exists in database, skip those:');
    disp(en(i)); % display exp
    pExp(i) = []; % exclude exp
    if isempty(pExp)
        disp('no new experiments to add');
        return
    end
end
%--------------------------------------------------------------------------
% group names       
%--------------------------------------------------------------------------
[~,pG] = cellfun(@dircontent,pExp,'UniformOutput',0);
pG = celltakeout(pG); % take out group names paths
[pG,groupnames] = takeout_tempfiles(pG); % take out temp files
a = regexpi(groupnames,'^[A-Z]{1,}\d{1,}'); % validate names
a(cellfun(@isempty,a)) = {0};
v = cell2mat(a);
numel(v)
if sum(v) == numel(groupnames) % if names are valid
    fprintf('Validated: group names\n');
else
    fprintf('Error: group names below are incorrect: \n');
    disp(char(groupnames(~v))); % display
    error('Please correct above group names');
end
%--------------------------------------------------------------------------
% MWT names 
%--------------------------------------------------------------------------
[~,pMWT] = cellfun(@dircontent,pG,'UniformOutput',0);
pMWT = celltakeout(pMWT); % get MWT path
[pMWT,~] = takeout_tempfiles(pMWT); % take out temp files
mwtname = cellfun(@dircontent,pG,'UniformOutput',0);
mwtname = celltakeout(mwtname); % get MWT path
%--------------------------------------------------------------------------
% unzip if zipped
%--------------------------------------------------------------------------
z = regexpcellout(mwtname,'[.]zip$');
pMWTzipped = pMWT(z);
fprintf('%d/%d MWT files zipped\n',sum(z),numel(pMWT));
if sum(z) > 0
   fprintf('unzipping:\n');
    for zi = 1:numel(pMWTzipped)
        processIntervalReporter(numel(pMWTzipped),1,'unzipping',zi)
        punizp = pMWTzipped{zi};
        p2unzip = fileparts(punizp);
        unzip(punizp,p2unzip);
        delete(punizp);
    end
    fprintf('unzip; done\n');
end
clear punizp p2unzip z
%**************************************************************************
% INTEGRATE TO MAIN DATABASE 
%**************************************************************************
% move new files to database 
%--------------------------------------------------------------------------
[~,pExp] = dircontent(pNewData);
for i = 1:numel(pExp)
    ps = pExp{i};
    pd = pDataBase;
    movefile(ps,pd,'f');
end
%--------------------------------------------------------------------------
% update database 
%--------------------------------------------------------------------------
MWTDB = makeMWTDatabase4(pDataBase,pDataBase);
%--------------------------------------------------------------------------
% reporting
%--------------------------------------------------------------------------
fprintf('integration of new data completed \n\n');
%--------------------------------------------------------------------------
end



%% SUBFUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [varargout] = dircontent(p,varargin)
% a = full dir, can be folder or files, b = path to all files, 
% c = only folders, d = paths to only folders

%% get dir
switch nargin
    case 0
        varargout{1} = [];
        return
    
    case 1
        % get dir
        [a,b,c,d] = dirc(p);
        varargout{1} = a;
        varargout{2} = b;
        varargout{3} = c;
        varargout{4} = d;
        return
    
    
    case 2
        ext = varargin{1};
        % dircontnetext
        cd(p); % go to directory
        a = {};
        a = dir(ext); % list content
        a = {a.name}'; % extract folder names only
        a(ismember(a,{'.','..','.DS_Store'})) = []; 
        b = {};
        for x = 1:size(a,1); % for all files 
            b{x,1} = [p,'/',a{x,1}]; % make path for files
        end
         varargout{1} = a;
         varargout{2} = b;
         return
    
         
    case 3
        if strcmp(varargin{1},'Option')
            option = varargin{2};
        end
        
end


switch option
    case 'MWT'
        p1 = regexp(genpath(p),':','split')';
        [~,fn2] = cellfun(@fileparts,p1,'UniformOutput',0);
        i = ismember(fn2,{''});
        fn2(i) = []; 
        p1(i) = [];
        i = regexpcellout(fn2,'\<\d{8}[_]\d{6}\>');
        varargout{1} = fn2(i);
        varargout{2} = p1(i);
    case 'MWTall'
        p1 = regexp(genpath(p),':','split')';
        [~,fn2] = cellfun(@fileparts,p1,'UniformOutput',0);
        i = ismember(fn2,{''});
        fn2(i) = []; 
        p1(i) = [];
        i = regexpcellout(fn2,'\<\d{8}[_]\d{6}');
        varargout{1} = fn2(i);
        varargout{2} = p1(i);

    otherwise
        error 'No such option';
end
end

function [a,b,c,d] = dirc(p)
% a = full dir, can be folder or files, b = path to all files, 
% c = only folders, d = paths to only folders
cd(p); % go to directory
a = {}; % create cell array for output
a = dir; % list content
a = {a.name}'; % extract folder names only
a(ismember(a,{'.','..','.DS_Store'})) = []; 
b = {};
c = {};
d = {};
for x = 1:size(a,1) % for all files 
    b(x,1) = {strcat(p,'/',a{x,1})}; % make path for files
    if isdir(b{x,1}) ==1 % if a path is a folder
        c(end+1,1) = a(x,1); % record name in cell array b
        d(end+1,1) = b(x,1); % create paths for folders
    else
    end
end
end

function [p2,names] = takeout_tempfiles(p)
    [~,names] = cellfun(@fileparts,p,'UniformOutput',0);
    i = cellfun(@isempty,names) | regexpcellout(names,'^[.]');
    names = names(~i);
    p2 = p(~i);        
    cellfun(@delete,p(i)); % delete temp file
end

function [A] = regexpcellout(C,searchterm,varargin)
    reglist = {'split','match'};
    % make sense of inputs
    switch nargin
        case 2
            option = 'logical';
            % validate first input as cell
            %if iscell(varargin{1}), C = varargin{1}; end
            %if ischar(varargin{2}), searchterm = varargin{2}; end
            B = regexp(C,searchterm);

        case 3
            %if iscell(varargin{1}), C = varargin{1}; end
            %if ischar(varargin{2}), searchterm = varargin{2}; end
            a = varargin{1};
            if ischar(a); option = a; end
            i = strcmp(reglist,option);
            if sum(i)==1
                B = regexp(C,searchterm,option);
            else
                B = regexp(C,searchterm);
            end

        otherwise
            error 'incorect number of inputs';
    end


    switch option
        case 'split'
            A = {};
            if iscell(B{1})
                for x = 1:numel(B)
                    col = size(B{x},2); 
                    A(x,1:col) = B{x}; 
                end
            else
                for x = 1:numel(B)
                    col = size(B(x),2); 
                    A(x,1:col) = B(x); 
                end
            end


        case'match'
            col = cell2mat(cellfun(@size,B,'UniformOutput',0));
            if size(col,2) > 1
                n = max(col(:,2));
            else
                n = 1;
            end
            A = cell(numel(B),n);
            for x = 1:numel(B)
                if isempty(B{x})==0
                    col = size(B{x,1},2);
                    A(x,1:col) = B{x};
                else
                    A(x,1) = {''}; 
                end
            end

        case 'logical'
             A = [];
            for x = 1:numel(B)
                if isempty(B{x})==0 
                    A(x,1) = B{x,1};
                else
                    A(x,1) = 0;
                end
            end
            A = logical(A); 
      
        otherwise
    end

end

function [A] = celltakeout(B,varargin)
    % if no option selection, examine input B
    switch nargin
        case 1

           % 2 level cell
            if iscell(B{1}) ==1
                i = cell2mat(cellfun(@size,B,'UniformOutput',0));
                onecolume = sum(i(:,2)==1) == numel(B); % see if only one colume
                multirow = sum(i(:,1)>1) ~= 0; % see if multiple rows
                if onecolume == 1 %onecolume == 1 && multirow == 1
                    option = 'multirow';
                elseif onecolume == 0 && multirow == 0
                    option = 'split';
                end
            end

            % 1 level cell
            if iscell(B{1}) == 0 && iscell(B) ==1
                if size(B,2) >1
                    % see if only one colume
                    onecolume = 0;
                else 
                    onecolume = 1;
                end
                if size(B,1) > 1
                    % see if multiple rows
                    multirow = 1;
                else 
                    multirow = 0;

                end

                if onecolume == 1 && multirow == 1
                    option = 'multirow';
                end

                if onecolume == 0 && multirow == 0
                    option = 'split';
                end
            end

        case 2
            option = varargin{1};

    end

    % if option is given
    switch option
        case 'onecolmultirow'

        case 'split'
            A = {};
            for x = 1:numel(B)
                col = size(B{x},2); 
                A(x,1:col) = B{x}; 
            end
        case 'multirow'
            A = {}; 
            for x = 1:numel(B)
                A = [A;B{x}]; 
            end
        case 'singlerow'
            A = {};
            for x = 1:numel(B)
                if isempty(B{x})==0; A(x,1) = B{x,1};
                else A(x,1) = {''}; end
            end
        case'match'
            A = {};
            for x = 1:numel(B)
                if isempty(B{x})==0; A(x,1) = B{x,1};
                else A(x,1) = {''}; end
            end
        case 'singlenumber'
            A = [];
            for x = 1:numel(B)
                if isempty(B{x})==0; A(x,1) = B{x,1};
                else A(x,1) = 0; end
            end
        case 'logical'
             A = [];
            for x = 1:numel(B)
                if isempty(B{x})==0; A(x,1) = B{x,1};
                else A(x,1) = 0; end
            end
            A = logical(A);

        otherwise
            error('no %s function',option);
    end

end

function processIntervalReporter(nfiles,itrgap,name,it_current)

    if ismember(it_current,1:itrgap:nfiles) == 1
        fprintf('%s: %d/%d\n',name,it_current,nfiles);
    end
end

function MWTDB = makeMWTDatabase4(pData,pMWTDB)
%% makeMWTDatabase(pData,pMWTDB)
% INPUT:
%     pData = path to MWT data
%         pData = '/Volumes/COBOLT/MWT';
%     pMWTDB = path to MWT database storage site
%         pMWTDB = '/Volumes/COBOLT/MWT';

% run these scripts %%
% pData = '/Volumes/COBOLT/MWT';
% pMWTDB = '/Volumes/COBOLT/MWT';
% MWTDB = makeMWTDatabase3(pData,pMWTDB);


%% path
if nargin <2
   pMWTDB = pData; 
end

%% function path
addpath('/Users/connylin/Dropbox/Code/Matlab/Library/General');


%% load current MWTDB
cd(pMWTDB);
MWTDB_current = load('MWTDB.mat','MWTDB');
MWTDB_old = MWTDB_current;

%% archive old MWTDB
% make new name
archivename = sprintf('MWTDB_%s',generatetimestamp);
% save old file as new name
eval(sprintf('%s = MWTDB_current.MWTDB;',archivename))
% archive
cd(pMWTDB);
eval(sprintf('save(''MWTDB.mat'',''%s'',''-append'')',archivename))


%% survey database
fprintf('\nGetting info from database %s... ',pData);
pMWTD = getMWTfile(pData);
fprintf('done\n');


%% GENERATE DATABASE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
D = struct;
% MWT index
T = parseMWTinfo(pMWTD);
D.text = T;


%% group name
fprintf('getting groupname... ');
% unique group names
fGU = unique(D.text.groupname);
T= table;
T.groupname_id = [1:numel(fGU)]';
T.groupname = fGU;
T.expN = nan(numel(fGU),1);
T.mwtN = nan(numel(fGU),1);
for x = 1:numel(fGU)
    i = ismember(D.text.groupname,fGU{x}); 
    a = unique(D.text.expname(i));
    T.expN(x) = numel(a);
    T.mwtN(x) = sum(i);
end
T.problem = false(numel(fGU),1);
T.problem(~regexpcellout(fGU,'^[A-Z]{1,}\d{1,}') & ~regexpcellout(fGU,'\<Test\>')) = true;
D.groupname = T;
% report problem
if sum(T.problem == true) > 1
    fprintf('\nproblem group names:\n');
    fGproblem = T.groupname(T.problem == true);
    disp(char(fGproblem));
    fprintf('\nfrom exp:\n');
    disp(char(unique(D.text.expname(ismember(D.text.groupname,fGproblem)))))
else
    fprintf('\ngroup names ok\n');
end


%% strain
fprintf('getting strain... ');
strain = D.text.strain;
strainU = unique(strain);
% validate
i = cellfun(@isempty,D.text.strain);
if sum(i) > 1
    fprintf('\nproblem strain names in experiment:\n');
    disp(char(unique(D.text.expname(i))));
    % display problem group names
    fprintf('\nproblem group names:\n');
    disp(char(unique(D.text.groupname(i))))
    % eliminate not qualified
    i = ~regexpcellout(strainU,'\<[A-Z]{1,}\d{1,}\>');
    strainU(i) = [];
else
    fprintf('\nstrain names ok\n');
end
fprintf('total %d valid strains\n',numel(strainU));
T = table;
T.strain_id = [1:numel(strainU)]';
T.strain = strainU;
D.strain = T;


%% rx
b = regexpcellout(D.text.rx,'_','split');
rx = {};
for x = 1:numel(b)
   rx = [rx;b{x}];
end
rxU = unique(rx);
T = table;
T.rxparts_id = [1:numel(rxU)]';
T.rxparts = rxU;
T.problem = false(numel(rxU),1);
rxProblemIgore = {'0A','0B','2A','2B','3A','3B'};
i = ismember(rxU,rxProblemIgore);
T.problem(i) = true;

T.expN = nan(numel(rxU),1);
T.mwtN = nan(numel(rxU),1);
T.solution = cell(numel(rxU),1);
for x = 1:numel(rxU)
    i = ~cellfun(@isempty,regexp(D.text.groupname,['_',rxU{x}],'once'));
    a = unique(D.text.expname(i));
    T.expN(x) = numel(a);
    b = D.text.groupname(i);
    T.mwtN(x) = numel(b);
end
D.rxparts = T;


%% rc
fGU = unique(D.text.rc);
T = table;
T.rc_id = [1:numel(fGU)]';
T.rc = fGU;
T.expN = nan(numel(fGU),1);
T.mwtN = nan(numel(fGU),1);
for x = 1:numel(fGU)
    i = ismember(D.text.rc,fGU{x}); 
    a = unique(D.text.expname(i));
    T.expN(x) = numel(a);
    T.mwtN(x) = sum(i);
end
D.rc = T;


%% cond
rxU = unique(D.text.rx);
rxU = rxU(2:end);
T = table;
T.rx_id = [1:numel(rxU)]';
T.rx = rxU;
T.expN = nan(numel(rxU),1);
T.mwtN = nan(numel(rxU),1);
for x = 1:numel(rxU)
    i = ismember(D.text.rx,rxU{x}); 
    a = unique(D.text.expname(i));
    T.expN(x) = numel(a);
    T.mwtN(x) = sum(i);
end
D.rx = T;


%% exp report
a = cellfun(@strcat,D.text.rc, cellfunexpr(D.text.rc,'_'),D.text.groupname,'UniformOutput',0);
A = unique(a);
gn = regexpcellout(A,'(?<=[_])\w{1,}','match');
r = regexpcellout(A,'_','split');
r = r(:,1);
T = table;
T.rc_groupname_id = [1:numel(A)]';
T.rc_groupname = A;
T.rc = r;
T.groupname = gn;
T.expN = nan(numel(A),1);
T.mwtN = nan(numel(A),1);
for x = 1:numel(A)
    i = ismember(D.text.rc,r{x}) & ismember(D.text.groupname,gn{x}); 
    a = unique(D.text.expname(i));
    T.expN(x) = numel(a);
    T.mwtN(x) = sum(i);
end
D.rc_groupname = T;


%% strain list - genotype +++++++++++++++++++++++++++++++++++++++++++++++++
% update strain/genotype list --------------------------------------------
SL = MWTDB_current.MWTDB.strain; % get old table
A = outerjoin(SL,D.strain,'Keys',{'strain'},'MergeKeys',1); % merge table
A.Properties.VariableNames{1} = 'strain_id'; % change id name
A.strain_id_right = [];
A = sortrows(A,'strain_id'); % sort table by strain id
id = find(isnan(A.strain_id)); % find new strains without id
A.strain_id(id) = id; % add new strain id 
% obtain new genotypes
fprintf('\n\n');
strain_noID = A.strain(id);
for x = 1:numel(strain_noID)
    strain_name = strain_noID{x}; % get strain name
    s = sprintf('enter genotype for %s: ',strain_name); % create prompt string
    a = input(s,'s'); % ask for genotype
    A.genotype{id(x)} = a; % enter genotype to database
end
D.strain = A;
%--------------------------------------------------------------------------
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

%% genotype
% find genes from genotype
a = D.strain.genotype;
b = regexp(a,'([a-z]+(-)\d+)|(wildtype*)|(?)','match');
% deal with exception (manual fix)
exception = {'C18H7.2'};
for x = 1:numel(exception)
    b(regexpcellout(a,exception{x})) = {exception(x)};
end
b(cellfun(@isempty,b)) = [];
% visulize what's not coded
ng = a(cellfun(@isempty,b)); 
disp(ng);
% get list of genes
b = celltakeout(b);
b = reshape(b,numel(b),1);
b(cellfun(@isempty,b)) = [];
b = unique(b);
b = [b;ng];
% merge new and old database
NEW = table(b,'VariableNames',{'gene'}); % get new genes
OLD = MWTDB_current.MWTDB.gene; % get old table
A = outerjoin(OLD,NEW,'Keys',{'gene'},'MergeKeys',1); % merge table
A = sortrows(A,'gene_id'); % sort table by gene id
id = find(isnan(A.gene_id)); % find new gene without id
A.gene_id(id) = id; % add new gene id 
D.gene = A; % add to main database


%% allele
a = D.strain.genotype; % get genotype
b = regexp(a,'([a-z]+(-)\d+)|(wildtype*)|(?)','match'); % get allele
% deal with exception
exception = {'C18H7.2'};
for x = 1:numel(exception)
    b(regexpcellout(a,exception{x})) = {exception(x)};
end
allele = cell(size(b));
for x = 1:size(b,1)
    c = b{x};
    allele{x} = cell(size(c));
    for y = 1:numel(c)
       if strcmp(c{y},'wildtype') == 1
           str = sprintf('(?<=(%s)\\()\\w+(?=\\))',c{y});
       else
          str = sprintf('(?<=(%s)\\()\\w+\\d+(?=\\))',c{y});
       end
       d = regexpcellout(a(x),str,'match');       
       allele{x}(y) = d;
    end
end
% create unique new allele list
b = allele;
b = celltakeout(b);
b = reshape(b,numel(b),1);
b(cellfun(@isempty,b)) = [];
b = unique(b);
b = [b;ng];
% merge new and old database
OLD = MWTDB_current.MWTDB.allele; % get old table
NEW = table(b,'VariableNames',{'allele'}); % get new genes
A = outerjoin(OLD,NEW,'Keys',{'allele'},'MergeKeys',1); % merge table
A = sortrows(A,'allele_id'); % sort table by gene id
id = find(isnan(A.allele_id)); % find new gene without id
A.allele_id(id) = id; % add new gene id 
D.allele = A; % add to main database


%% allele - gene pair (suspend)
% OLD = MWTDB_current.MWTDB.gene_allele_id; % get old table
% 
% %%
% genes = b;
% genes = celltakeout(genes);
% a = cell(numel(genes),2);
% for x = 1:numel(genes)
%     a(x,1) = genes(x);
%     a(x,2) = allele(x);
% end
% a(cellfun(@isempty,a(:,2)),:) = [];
% % translate to id
% genes = a(:,1);
% allele = a(:,2);
% i = ismember(genes,'?');
% [i,j] = ismember(allele,D.allele.allele);
% allele_id = j;
% [i,j] = ismember(genes,D.gene.gene);
% gene_id= j;
% a = [gene_id allele_id];
% a = unique(a,'rows');
% t = table;
% t.gene_id = gene_id;
% t.allele_id = allele_id;
% D.gene_allele_id = t;


%% gene sets
% exp set based on genes
% %% archive
a= {...
    'G protein',{'acy-1';'dgk-1';'egl-30';'egl-8';'goa-1';'gsa-1';'kin-2';'pde-4';'pkc-1';'ric-8'};
    'calcium',{'cnb-1'};
    'dopamine',{'cat-2', 'dop-1'};
    'synaptic',{'eat-16';'tom-1';'unc-13'};
    'nicotinic ion channel',{'eat-2'};
    'VG K channel',{'egl-2';'unc-103'};
    'neuropeptide Y',{'npr-1';'flp-18';'flp-20';'flp-21';'wildtyle;hawaiian'};
    'Glutamate',{'glr-1'};
    'neuroligin',{'mir-1';'nlg-1';'nrx-1';'pros-1'};
    'etoh insensitive',{'ptr-6'};
    'BK',{'slo-1'}
    };
D.expset = cell2table(a,'VariableNames',{'expset','genes'});


%% save
cd(pMWTDB);
MWTDB = D;
save('MWTDB.mat','MWTDB');
% writetable(D.rc_groupname,'rc_groupname.csv')
% writetable(D.text,'mwt.csv'); 


end

function [pMWT,fMWT,pMWTBadName] = getMWTfile(p,type)
%% getMWTfile
% [pMWT,fMWT] = getMWTfile(p)
% input
%     p = path of folder contining experiment folders
%     type = 
%         regular - go through regular process
%         zip - the p folder contains zipped mwt files, find zip files
% output
%     pMWT = path to MWT zip files or folders
%     fMWT = names of MWT zip files or folder
%     pMWTBadName = mwt with bad names

%% process input
if nargin <2
   type = 'regular'; 
end

%% get fMWT and pMWT
switch type
    case 'regular'
        [~,~,~,pEf] = dircontent(p); % get exp path
        % get group folder
        [~,~,~,pG] = cellfun(@dircontent,pEf,'UniformOutput',0);
        pG = celltakeout(pG);
        % get mwt paths
        [fMWT,pMWT] = cellfun(@dircontent,pG,'UniformOutput',0);
        pMWT = celltakeout(pMWT,'multirow');
        fMWT = celltakeout(fMWT,'multirow');
    case 'zip'
        [pMWT,fMWT] = getMWTzipped(p);
end

%% process output
% flag problem mwt folder names
j = regexpcellout(fMWT,'\<\d{8}[_]\d{6}'); % find files with MWT file name as a start
% validate MWT folders
i = regexpcellout(fMWT,'(\<\d{8}[_]\d{6}\>)|(\<\d{8}[_]\d{6}[.]zip\>)'); % find files with proper MWT names
if sum(j) ~= sum(i)
    pMWTBadName = pMWT(j == true & i == false);
end

pMWT(~i) = [];
fMWT(~i) = [];

end

function T = parseMWTinfo(pMWTD)

%% accomodate for empty entry
if isempty(pMWTD) == true
    T = [];
else
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
    % T.exp_date = cellfun(@str2num,regexpcellout(fE,'\<\d{8}','match'));
    T.exp_date = cellfun(@str2num,regexpcellout(fE,'^\d{8}','match'));

    T.tracker = regexpcellout(fE,'(?<=\<\d{8})[A-Z]','match');
    T.expter = regexpcellout(fE,'(?<=\<\d{8}[A-Z][_])[A-Z]{2}','match');
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
    
end
end







