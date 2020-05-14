
%% input variables
pExpFolder = '/Users/connylin/Dropbox/Lab/MWT_Analysis_ByExp/Recovery 10sISI';
pAnalysis = '/Users/connylin/Dropbox/Lab/MWT_Analysis';
pData = '/Volumes/ParahippocampalGyrus/MWT/Data';
% add function
addpath('/Users/connylin/Dropbox/MATLAB/Function_Library_Public');
% create experiment set folder
% pDest = ['/Users/connylin/Dropbox/Lab/MWT_Analysis_ByExp/',expsetName];
% if isdir(pDest) == 0; mkdir(pDest); end


%% get group info
% rasterData: row = worm, column = response in every .2seconds
fG = 'N2_400mM';
fD = 'rasterPlot_95_105_N_471.rasterData';
p = [pExpFolder,'/',fG,'/',fD];
d = dlmread(p);
% response around the tap
t1 = 25;
t2 = 32;
% find ones has no reversals
i = sum(d(:,t1:t2),2);
x = 1:numel(i);
% y = zeros(numel(i),1);
% y(i) = 10;
y = i;
plot(x,y)


%%
plateInfo = readtable([pExpFolder,'/plate_info.csv']);
% reconstruct path
pMWT = cell(size(plateInfo,1),1);
for x = 1%:size(plateInfo,1)
    % get trinitySummary.mat file
    p = sprintf('%s/%s/%s/%s/%s',pAnalysis,plateInfo.expname{x},plateInfo.groupname{x},plateInfo.mwtname{x},'trinitySummary.mat');
    load(p)
end



% 
% expsetName = 'N2 45ISI';
% expSearchTerm = '100s30x45s10s';
% groupSearchTerm = '(\<N2_400mM\>|\<N2\>)';
% % search type:
% %     1 = only groupSearchTerm and expSearchTerm
% %     2 = not group mentioned and within exp name specified
% %     3 = all groups within experiments that contains the group terms
% searchType = 1;
% % time
% tStart = 85;
% ISI = 45;
% tapN = 30;
% tEnd = 100+(ISI*tapN)+10;
% tInt = 10;
% % calculate time inputs
% t = [tStart tStart+tInt:ISI:tEnd];
% timePairs = [t;t+tInt]';
% disp('time pairs:');
% disp(timePairs);


%% Paths




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
        fprintf('\nGroups included [%d]:\n',numel(unique(fG)));
        disp(char(unique(fG)));
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
pMWTC = pMWT(cellfun(@numel,cellfun(@dircontent,pMWTA,cellfunexpr(pMWTA,'trinitySummary.mat'),'UniformOutput',0)) == 0);
fprintf('MWT to chor: %d ',numel(pMWTC));
% confirm
if input('Stop [1] or Press any other key to continue: ') == 1 ;  return; end

% chor
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
pMWTA(ismember(pMWTA,pMWTCA)) = [];
pMWTA = [pMWTA;pMWTCAV];
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
        [figure1,plateSumm,savename] = rasterPlot_colorSpeed(pathGrouped{gi},timePairs(ti,1),timePairs(ti,2),'NWORMS',Inf,'InputType',2);
        % create save folder
        pSave = [pDest,'/',fGU{gi},'/rasterPlot_colorSpeed'];
        if isdir(pSave) == 0; mkdir(pSave); end
        cd(pSave);
        set(figure1,'PaperPositionMode','auto'); % set to save as appeared on screen
        print (figure1,'-depsc', '-r1200', savename); % save as eps
        save([savename,'.rasterData'],'plateSumm','-ascii');
        close;
    end
end


%% report done
fprintf('* %s  done *\n',mfilename)
















