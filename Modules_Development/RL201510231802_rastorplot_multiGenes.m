function RL201510231802_rastorplot_multiGenes(pMWT,expsetName,varargin)
%% input variables
% expsetName = 'slo1 10sISI';
% searchTerm = '100s30x10s10s';
% groupSearchTerm = '(NM1968)|(NM1630)|(BZ142)';
%     OPTION:
%         'checkStop':
%             0 = default, go through the whole program
%             1 = check groups to be analyzed
% 

%% default setting
tStart = 85;
tEnd = 405;
tInt = 10;
checkStop = 0;

%% deal with varagin
if nargin > 3
    A = varargin;
    if isinteger(numel(A)/2) == 1
        error('variable entries incorrect')
    end
    callName = A([1:2:numel(A)]);
    setValue = cell2mat(A([2:2:numel(A)]));
    for x = 1:numel(callName)
        str = sprintf('%s = %d;',callName{x},setValue(x));
        eval(str)
    end
end

%% deal with setting variables
% calculate time inputs
timePairs = [tStart:tInt:tEnd-tInt; tStart+tInt:tInt:tEnd]';



%% Paths
pAnalysis = '/Users/connylin/Dropbox/Lab/MWT_Analysis';
pData = '/Volumes/ParahippocampalGyrus/MWT/Data';
% add function
addpath('/Users/connylin/Dropbox/MATLAB/Function_Library_Public');
% create experiment set folder
pDest = ['/Users/connylin/Dropbox/Lab/MWT_Analysis_ByExp/',expsetName];
if isdir(pDest) == 0; mkdir(pDest); end




%% chor Trinity
addpath('/Users/connylin/Dropbox/MATLAB/programs_rankinlab/Library/Modules/Chor');
[pMWTA,fval,chorscript] = chormaster3('Trinity', pMWT,pAnalysis,pData,'old');



%% get path to experiments
% check if all has analysis output
fprintf('\n\nvalidating trinity output...\n');
fnameVal = fval{1};
[fn,p] = cellfun(@dircontent,pMWTA,cellfunexpr(pMWTA,fnameVal),'UniformOutput',0);
nFile = cellfun(@numel,fn);
nVal = nFile > 0;
if sum(nVal) ~= numel(pMWTA)
    warning('not all MWT has %s output',fnameVal);
    % eliminate ones without output and prepare report
    [pG,fMWT] = cellfun(@fileparts,pMWTA,'UniformOutput',0);
    fprintf('the following MWT plates are removed from analysis:\n');
    disp(char(fMWT(~nVal)))
    pMWTA = pMWTA(nVal);
else
    fprintf('all %d MWT files have %s outputs\n',numel(pMWTA),fnameVal);
end

% find mwt file info
[pG, fMWT] = cellfun(@fileparts,pMWTA,'UniformOutput',0);
[pExp, fG] = cellfun(@fileparts,pG,'UniformOutput',0);
[~, fE] = cellfun(@fileparts,pExp,'UniformOutput',0);



%% create outputs
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



%% [Disable] copy to sort into groups
% display('Copying files');
% fGU = unique(fG);
% for fGUi = 1:numel(fGU)
%     % get folder name
%     pDestG = [pDest,'/',fGU{fGUi}];
%     fprintf('group %d/%d: %s\n',fGUi,numel(fGU),fGU{fGUi});
%     if isdir(pDestG) == 0; mkdir(pDestG); end
%     i = ismember(fG,fGU(fGUi));
%     pT = pMWTA(i);
%     for pTi = 1:numel(pT)
%         ps = pT{pTi};
%         [~,fn] = fileparts(ps);
%         pd = [pDestG,'/',fn];
%         if isdir(pd) == 0; mkdir(pd); end
%         [~,psfiles] = dircontent(ps);
%         pdfiles = regexprep(psfiles,ps,pd);
%         if numel(pdfiles) > 70
%             cellfun(@copyfile,psfiles,cellfunexpr(psfiles,pd));
%         else
%             copyfile(ps,pd,'f') % this function won't work if too many files
%         end
%     end
% end



%% run raster plot (use pMWT)
addpath('/Users/connylin/Dropbox/MATLAB/Programs_RankinLab/Library/Modules/rasterPlot_colorSpeed');


% group MWT paths into groups
[pG,fMWT] = cellfun(@fileparts,pMWTA,'UniformOutput',0);
[pE,fG] = cellfun(@fileparts,pG,'UniformOutput',0);
fGU = unique(fG);
pathGrouped = cell(size(fGU));
for x = 1:numel(fGU)
    pathGrouped{x} = pMWTA(ismember(fG,fGU{x}));
end

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
















