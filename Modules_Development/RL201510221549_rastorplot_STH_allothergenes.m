% add function
addpath('/Users/connylin/Dropbox/MATLAB/Function_Library_Public');
pAnalysis = '/Users/connylin/Dropbox/Lab/MWT_Analysis';
pData = '/Volumes/ParahippocampalGyrus/MWT/Data';
expsetName = 'slo1 10sISI';
searchTerm = '100s30x10s10s';
groupSearchTerm = '(NM1968)|(NM1630)|(BZ142)';
timePairs = [85:10:395; 95:10:405]';

%% search 
% search exp
[~,~,fE,pE] = dircontent(pData);
i = regexpcellout(fE,searchTerm);
fprintf('\nexperiments found [%d]:\n',sum(i));
disp(char(unique(fE(i))))
pE = pE(i);
fE = fE(i);

% search groups
[fG,pG] = cellfun(@dircontent,pE,'UniformOutput',0);
pG = celltakeout(pG);
fG = celltakeout(fG);
i = regexpcellout(fG,groupSearchTerm);
pG = pG(i);
fG = fG(i);
fprintf('\ngroups found [%d total, %d unique]:\n',sum(i), numel(unique(fG)));
disp(char(unique(fG)))

% get experiment with group names
[pE,fG] = cellfun(@fileparts,pG,'UniformOutput',0);
[~,fE] = cellfun(@fileparts,pE,'UniformOutput',0);
fprintf('\nExperiments included [%d]:\n',numel(unique(fE))); 
disp(char(unique(fE)));
pE = unique(pE);
[~,~,fG,pG] = cellfun(@dircontent,pE,'UniformOutput',0);
pG = celltakeout(pG); fG = celltakeout(fG);
fprintf('\nGroups included [%d]:\n',numel(unique(fG)));
disp(char(unique(fG)));


% get MWT files within criteria
[~,pMWT] = cellfun(@dircontent,pG,cellfunexpr(pG,'*.zip'),'UniformOutput',0);
pMWT = celltakeout(pMWT);

% reconfirm groups
[pG,fMWT] = cellfun(@fileparts,pMWT,'UniformOutput',0);
[pE,fG] = cellfun(@fileparts,pG,'UniformOutput',0);
[~,fE] = cellfun(@fileparts,pE,'UniformOutput',0);

% make sure MWT files are unique
if numel(unique(pMWT)) ~= numel(pMWT)
   error('mwt paths not unique'); 
end


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

% create experiment set folder
pDest = ['/Users/connylin/Dropbox/Lab/MWT_Analysis_ByExp/',expsetName];
if isdir(pDest) == 0; mkdir(pDest); end

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



%% copy to sort into groups
fGU = unique(fG);
for fGUi = 2:numel(fGU)
    % get folder name
    pDestG = [pDest,'/',fGU{fGUi}];
    if isdir(pDestG) == 0; mkdir(pDestG); end
    i = ismember(fG,fGU(fGUi));
    pT = pMWTA(i);
    for pTi = 19:numel(pT)
        ps = pT{pTi};
        [~,fn] = fileparts(ps);
        pd = [pDestG,'/',fn];
        if isdir(pd) == 0; mkdir(pd); end
        [~,psfiles] = dircontent(ps);
        pdfiles = regexprep(psfiles,ps,pd);
        cellfun(@copyfile,psfiles,cellfunexpr(psfiles,pd));
%         copyfile(ps,pd,'f') % this function won't work if too many files
    end
end



%% run raster plot
addpath('/Users/connylin/Dropbox/MATLAB/Programs_RankinLab/Library/Modules/rasterPlot_colorSpeed');
pE = ['/Users/connylin/Dropbox/Lab/MWT_Analysis_ByExp/',expsetName];
[~,~,fG,pG] = dircontent(pE);

for gi = 1:numel(pG)
    fprintf('Raster plotting group [%s]...\n',fG{gi});
    for ti = 1:size(timePairs,1)
        str = sprintf('rasterPlot_%d_%d*',timePairs(ti,1),timePairs(ti,2));
        fprintf('--%s',str);
        if numel(dircontent(pG{gi},str)) == 0
            plateSumm = rasterPlot_colorSpeed(pG{gi},timePairs(ti,1),timePairs(ti,2),'NWORMS',Inf);
        end
    end
    
end



%% report done
fprintf('* %s  done *\n',mfilename)
















