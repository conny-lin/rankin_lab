function [MWTzip, copy_error_dump] = consolidate_MWTFiles(pTarget,pKeep)
%% view current search path
addpath('/Users/connylin/Dropbox/Code/Matlab/Library/General');
addpath('/Users/connylin/Dropbox/rl/Code/Modules/MWTDatabase');


% pTarget = '/Volumes/IRONMAN/MWT data/Checked_Data_AnalysisOutput_Archive';
% pDatabase = '/Volumes/COBOLT/MWT';


% load database
load('/Users/connylin/Dropbox/rl/MWTDB/MWTDB.mat');
mwtname_database = MWTDB.text.mwtname;
mwtpath_database = MWTDB.text.mwtpath;

copy_error_dump = {};


%% get mwt folder path in pTarget
p = genpath(pTarget);
pSubFolders = regexp(p,':','split');
pSubFolders(cellfun(@isempty,pSubFolders)) = [];
% look into folders for MWT folders
[pMH,foldernames] = cellfun(@fileparts,pSubFolders,'UniformOutput',0);
i = regexpcellout(foldernames','(\<\d{8}[_]\d{6}\>)|');
pMWTfolders = pSubFolders(i);
pFolders = pSubFolders(~i);
mwtname_target = foldernames(i);
fprintf('%d MWT folder found\n',numel(mwtname_target));

%% get mwt folder path in pKeep
p = genpath(pKeep);
pSubFolders = regexp(p,':','split');
pSubFolders(cellfun(@isempty,pSubFolders)) = [];
% look into folders for MWT folders
[~,foldernames] = cellfun(@fileparts,pSubFolders,'UniformOutput',0);
i = regexpcellout(foldernames','(\<\d{8}[_]\d{6}\>)|');
pMWTfolders_keep = pSubFolders(i);
mwtname_keep = foldernames(i);
fprintf('%d MWT folder found in KEEP folder\n',numel(mwtname_keep));

%% check if MWT folder already exist in COBOLT
% check folder name
[j,k] = setdiff(mwtname_target, mwtname_database);
mwtname_target_Notdatabase = mwtname_target(k);

[j1,k1] = setdiff(mwtname_target_Notdatabase, mwtname_keep);
if isempty(k) == 1 || isempty(k1) == 1
    fprintf('no new MWT folder found\n');
else
    fprintf('move %d folder content to keep folder\n',numel(k1));
    psourcefolder = pMWTfolders(k1);
    pdestfolder = regexprep(psourcefolder,pTarget,pKeep);
    cellfun(@movefile,psourcefolder,pdestfolder);
end

[j,k] = ismember(mwtname_target, mwtname_database);
mwtpath_db_match = mwtpath_database(k(j));

%% check if new files within each folder
fprintf('copy content for each %d folders:\n',numel(mwtpath_db_match));

for x = 1:numel(pMWTfolders)
    ptg = pMWTfolders{x};
    [~,ftg] = fileparts(ptg);
    pdb = mwtpath_db_match{x};
    [~,fdb] = fileparts(pdb);
    if strcmp(ftg,fdb) ~= 1; error('mwt name match incorrect'); end
    
    [fpdb,~] = dircontent(pdb);
    [fptg, pptg] = dircontent(ptg);
    fi = ismember(fptg,fpdb);
    if sum(fi) ~= numel(fi)
        fprintf('mwt %d: %d files not found in database ',x,sum(~fi))
        % copy file to database
        psource = pptg(~fi);
        % exclude files
        str = '([.]swanlake2all[.])|([.]swanlake2.dat)|([.]swn)';
        i = regexpcellout(psource,str);
        psource(i) = [];
        fprintf(', %d files valid to tsf ',numel(psource))
        % get pdesignation
        pdest = regexprep(psource,ptg,pdb);

        % copy file
        [sc,msg,msgid] = cellfun(@copyfile,psource,pdest,'UniformOutput',0);
        sc = cell2mat(sc);
        if sum(sc) ~= sum(~fi)
            fprintf('mwt %d: %d files not copied correctly\n',x,sum(~sc));
            copy_error_dump = [copy_error_dump;psource(~sc)];
        else
            fprintf('-- copied\n');
        end
    end
end
fprintf('%d files did not copy correctly\n',numel(copy_error_dump));






%% look into none MWT folders to find zip files
str = '\<\d{8}[_]\d{6}[.]zip\>';
MWTzip = {};
for i = 1:numel(pFolders)
   p = pFolders{i};
   [~,pZip] = dircontent(p,'*.zip');
   if isempty(pZip) == 0
       k = regexpcellout(pZip,str);
       MWTzip = [MWTzip;pZip(k)];
   end
end
fprintf('%d zipeed mwt files found\n',numel(MWTzip));
%%
cellfun(@fileparts,MWTzip,'UniformOutput',0)

%% unzip all







