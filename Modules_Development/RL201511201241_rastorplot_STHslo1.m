clear;
%% input variables
expsetName = 'slo1 10sISI simple';
tStart = 85;
tEnd = 405;
tInt = 10;
timePairs = [tStart:tInt:tEnd-tInt; tStart+tInt:tInt:tEnd]';


%
%% Paths
pAnalysis = '/Volumes/COBOLT/MWT';
pData = '/Volumes/COBOLT/MWT';
% add function
addpath('/Users/connylin/Dropbox/MATLAB/Function_Library_Public');
% create experiment set folder
pDest = ['/Users/connylin/Dropbox/RL/MWT_Analysis_ByExp/',expsetName];
if isdir(pDest) == 0; mkdir(pDest); end


%% get database
addpath('/Users/connylin/Dropbox/RL/Code/Library/Modules/MWTDatabase');
MWTD = load(sprintf('%s/MWTDatabase.mat',pData));
MWTD = MWTD.MWTDatabase;
A = MWTD.mwt;

% find liquid transfer date
a = find(regexpcellout(A.groupname,'Liquid'));
liquidTsfDate = A.exp_date(a(end));

% find targets
i = regexpcellout(A.strain,'(NM1968)|(NM1630)|(BZ142)|(CX3940)') &...
    ismember(A.rx,{'NA','400mM'}) &...
    A.ISI == 10 &...
    A.preplate == 100 &...
    A.tapN == 30 &...
    A.exp_date > liquidTsfDate;
D = MWTD.mwt(i,:);

% find experiments with controls
expT = unique(D.expname);
i = ismember(A.expname,expT) & ...
    ismember(A.strain,{'N2'}) & ...
    ismember(A.rx,{'NA','400mM'}) &...
    A.ISI == 10 &...
    A.preplate == 100 &...
    A.tapN == 30 &...
    A.exp_date > liquidTsfDate;

D = [D;A(i,:)];

% eliminate duplicates
[~,i] = unique(D.mwt_id);
D(~i,:) = [];

% display names
display('group names:');
disp(char(unique(D.groupname)));


% create outputs
writetable(D,[pDest,'/plate_info.csv']);
% make table
save([pDest,'/plate_info.mat'],'D');



%% chor Trinity
% look if all files have required files
reqFile = 'trinitySummary.mat';
pMWT = D.mwtpath;
val = false(size(pMWT));
for mwti = 1:numel(pMWT)
    a = dircontent(pMWT{mwti},reqFile);
    if isempty(a) == 0; val(mwti) = true; end
end
fprintf('%d/%d MWT do not have %s file\n',sum(~val), numel(val),reqFile);

% check if no reqFile has summary file
reqFile = '*trinity*.dat';
pMWT = pMWT(~val);
val = false(size(pMWT));
for mwti = 1:numel(pMWT)
    a = dircontent(pMWT{mwti},reqFile);
    if isempty(a) == 0; val(mwti) = true; end
end
fprintf('%d/%d MWT have %s file\n',sum(val), numel(val),reqFile);

% chor
pMWT = pMWT(~val);
addpath('/Users/connylin/Dropbox/RL/Code/Library/Modules/Chor');
chormaster4('Trinity',pMWT);



%% run raster plot (use pMWT)
% look if all files have required files
reqFile = 'trinitySummary.mat';
pMWT = D.mwtpath;
val = false(size(pMWT));
for mwti = 1:numel(pMWT)
    a = dircontent(pMWT{mwti},'trinitySummary.mat');
    if isempty(a) == 0; 
        val(mwti) = true; 
    else
        a = dircontent(pMWT{mwti},'*trinity*.dat');
        if isempty(a) == 0
            val(mwti) = true; 
        end
    end
end
fprintf('%d/%d MWT do not have %s file\n',sum(~val), numel(val),reqFile);
D(~val,:) = [];
% create outputs
writetable(D,[pDest,'/plate_info.csv']);
% make table
save([pDest,'/plate_info.mat'],'D');
% group paths into groups
fG = D.groupname;
fGU = unique(D.groupname);
pathGrouped = cell(size(fGU));
for x = 1:numel(fGU)
    pathGrouped{x} = D.mwtpath(ismember(fG,fGU{x}));
end

% run
addpath('/Users/connylin/Dropbox/RL/Code/Library/Modules/rasterPlot_colorSpeed');

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


% report done
fprintf('* %s  done *\n',mfilename)
















