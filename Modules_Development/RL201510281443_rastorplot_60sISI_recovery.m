%% input variables
expsetName = 'Recovery 60sISI';

expSearchTerm = '100s31x60s10s';
groupSearchTerm = '(\<N2\>|\<N2_400mM\>)';
groupExclude = '(\<NM1630_400mMx3dx12h_400mMx4d\>|\<N2_400mMx3dx12h_400mMx4d\>|\<NM1630\>|\<NM1630_400mM\>)';
% search type:
%     1 = only groupSearchTerm and expSearchTerm
%     2 = not group mentioned and within exp name specified
%     3 = all groups within experiments that contains the group terms
%     4 = only experiments containing all groupSearchTerm, and exclude
%     groups specified in groupExclude
searchType = 4;
% paths
pAnalysis = '/Users/connylin/Dropbox/Lab/MWT_Analysis';
pData = '/Volumes/ParahippocampalGyrus/MWT/Data';
% time
tStart = 85;
ISI = 60;
tapN = 30;
% tEnd = 100+(ISI*tapN);
% tEnd = 695; % 10sISI recovery
tEnd = 100+(ISI*tapN)+(5*60-5); % 60sISI recovery
tInt = 10;
% calculate time inputs
t = [tStart tStart+tInt:ISI:tEnd];
timePairs = [t;t+tInt]';
disp('time pairs:');
disp(timePairs);


%% Paths
% add function
addpath('/Users/connylin/Dropbox/MATLAB/Function_Library_Public');
% create experiment set folder
pDest = ['/Users/connylin/Dropbox/Lab/MWT_Analysis_ByExp/',expsetName];
if isdir(pDest) == 0; mkdir(pDest); end


%% search 
[~,~,fE,pE] = dircontent(pData);
i = regexpcellout(fE,expSearchTerm);
fprintf('\nexperiments found with search term [%s]: %d\n',expSearchTerm,sum(i)); disp(char(unique(fE(i))));
pE = pE(i); fE = fE(i);

[~,~,fG,pG] = cellfun(@dircontent,pE,'UniformOutput',0);
pG = celltakeout(pG); fG = celltakeout(fG);
switch searchType
    case 1 % only groups specified
        i = regexpcellout(fG,groupSearchTerm);
        pG = pG(i); fG = fG(i);

    
    case 2 % groups exept for groups specified
        i = cellfun(@numel,regexp(fG,groupSearchTerm)) ==0;
        pG = pG(i); fG = fG(i);

    
    case 3 % all groups within experiments containing the group
        i = regexpcellout(fG,groupSearchTerm);
        pG = pG(i); fG = fG(i);
        % get experiments containing specified group names
        [pE,fG] = cellfun(@fileparts,pG,'UniformOutput',0);
        [~,fE] = cellfun(@fileparts,pE,'UniformOutput',0);
        fprintf('\nExperiments included [%d]:\n',numel(unique(fE))); 
        disp(char(unique(fE)));
        pE = unique(pE);
        [~,~,fG,pG] = cellfun(@dircontent,pE,'UniformOutput',0);
        pG = celltakeout(pG); fG = celltakeout(fG);
    
        
    case 4 % only experiments with all groups found
        i = regexpcellout(fG,groupSearchTerm);
        pG = pG(i); fG = fG(i);
        % get experiments containing specified group names
        [pE,fG] = cellfun(@fileparts,pG,'UniformOutput',0);
        [~,fE] = cellfun(@fileparts,pE,'UniformOutput',0);
        fprintf('\nExperiments included [%d]:\n',numel(unique(fE))); 
        disp(char(unique(fE)));
        pE = unique(pE);
        [~,~,fG,pG] = cellfun(@dircontent,pE,'UniformOutput',0);
        pG = celltakeout(pG); fG = celltakeout(fG);
        i = cellfun(@numel,regexp(fG,groupExclude)) ==0;
        pG = pG(i); fG = fG(i);
        
end
fprintf('\nGroups found: [%d total, %d unique]\n',sum(i), numel(unique(fG)));
disp(char(unique(fG)));


%% get MWT files within criteria
[~,pMWT] = cellfun(@dircontent,pG,cellfunexpr(pG,'*.zip'),'UniformOutput',0);
pMWT = celltakeout(pMWT);
% make sure MWT files are unique
if numel(unique(pMWT)) ~= numel(pMWT); error('mwt paths not unique'); end
% report number of MWT files
fprintf('\nMWT files included: %d\n\n',numel(pMWT));
% report MWT files per group
[~,fG] = cellfun(@fileparts,cellfun(@fileparts,pMWT,'UniformOutput',0),'UniformOutput',0);
tabulate(fG)


%% chor Trinity
% choose files that does not contain Trinity files
pMWTA = regexprep(regexprep(pMWT,pData,pAnalysis),'(.zip)','');
val = false(size(pMWTA)); 
% check if has analysis directory
for x = 1:numel(pMWTA)
    if isdir(pMWTA{x}) == false; 
        val(x) = true; 
    elseif numel(dircontent(pMWTA{x},'trinitySummary.mat')) == 0
        val(x) = true; 
    end
end
pMWTC = pMWT(val);
fprintf('MWT to chor: %d/%d\n',numel(pMWTC),numel(pMWTA));

% chor
% confirm
if input('Stop [1] or Press any other key to continue: ') == 1 ;  return; end
addpath('/Users/connylin/Dropbox/MATLAB/programs_rankinlab/Library/Modules/Chor');
[pMWTCA,fval,chorscript] = chormaster3('Trinity', pMWTC,pAnalysis,pData,'old');


%% get path to valid MWT
% check if all has analysis output
fprintf('\n\nvalidating trinity output...\n');
fnameVal = fval{1};
[fn,p] = cellfun(@dircontent,pMWTCA,cellfunexpr(pMWTCA,fnameVal),'UniformOutput',0);
nFile = cellfun(@numel,fn);
nVal = nFile > 0;
if sum(nVal) ~= numel(pMWTCA)
    warning('not all MWT has %s output',fnameVal);
    % eliminate ones without output and prepare report
    [~,fMWT] = cellfun(@fileparts,pMWTCA,'UniformOutput',0);
    fprintf('the following MWT plates are removed from analysis:\n');
    disp(char(fMWT(~nVal)))
    pMWTCAV = pMWTCA(nVal);
else
    fprintf('all %d MWT files have %s outputs\n',numel(pMWTCA),fnameVal);
end

pMWTA = unique([pMWTA;pMWTCAV]);
% report MWT files per group
[~,fG] = cellfun(@fileparts,cellfun(@fileparts,pMWTA,'UniformOutput',0),'UniformOutput',0);
fprintf('\nMWT plates going for raster plot: \n');
tabulate(fG)



%% create outputs
% find mwt file info
[pG, fMWT] = cellfun(@fileparts,pMWTA,'UniformOutput',0);
[pExp, fG] = cellfun(@fileparts,pG,'UniformOutput',0);
[~, fE] = cellfun(@fileparts,pExp,'UniformOutput',0);

% make table
T = table;
T.mwtname = fMWT;
T.groupname = fG;
T.expname = fE;
cd(pDest);
writetable(T,'plate_info.csv');

% make table
T.pMWT = pMWTA;
save([pDest,'/plate_info.mat'],'T');


%% group MWT paths into groups
[pG,fMWT] = cellfun(@fileparts,pMWTA,'UniformOutput',0);
[pE,fG] = cellfun(@fileparts,pG,'UniformOutput',0);
fGU = unique(fG);
pathGrouped = cell(size(fGU));
for x = 1:numel(fGU)
    pathGrouped{x} = pMWTA(ismember(fG,fGU{x}));
end


%% run raster plot (use pMWT)
addpath('/Users/connylin/Dropbox/MATLAB/Programs_RankinLab/Library/Modules/rasterPlot_colorSpeed');
for gi = 1:numel(fGU)
    fprintf('Raster plotting group [%s]...\n',fGU{gi});
    for ti = 1:size(timePairs,1)
        % reporting
        str = sprintf('rasterPlot_%d_%d*',timePairs(ti,1),timePairs(ti,2));
        fprintf('--raster plot [%d/%d]\n',ti,size(timePairs,1));
        % run raster plot
        [figure1,~,savename] = rasterPlot_colorSpeed(pathGrouped{gi},timePairs(ti,1),timePairs(ti,2),'NWORMS',Inf,'InputType',2);
        % create save folder
        pSave = [pDest,'/',fGU{gi},'/rasterPlot_colorSpeed'];
        if isdir(pSave) == 0; mkdir(pSave); end
        cd(pSave);
        set(figure1,'PaperPositionMode','auto'); % set to save as appeared on screen
        print (figure1,'-depsc', '-r1200', savename); % save as eps
        close;
    end
end


%% report done
fprintf('* %s  done *\n',mfilename)
















