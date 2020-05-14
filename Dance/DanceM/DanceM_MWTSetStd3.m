function [MWTSet,pSave] = DanceM_MWTSetStd3(pMWT,funpath,VarIn,pSave,varargin)
% DanceM_MWTSetStd2 - use DanceM_createIndex2 for flexibility

%% default
useoldfile = false;
timestamp = 'off';
outputsuffix = '';
copyfun = 'on';
vararginProcessor;


%% save folder
funpath = [funpath,'.m'];
[~,funame] = fileparts(funpath);
pSave = [pSave,'/',funame];
times = generatetimestamp;
if strcmp(timestamp,'on') 
   pSave = sprintf('%s_%s%s',pSave, outputsuffix, times);
end

% MWTSetval = false;
% if useoldfile
%     pMWTSet = sprintf('%s/%s/%s.mat',pSave,funame,funame);
%     if exist(pMWTSet)==2
%         display 'found prexisting mat file'
%         MWTSet = load(pMWTSet,'MWTSet'); MWTSet = MWTSet.MWTSet;
%         % check if pMWT matches
%         p1 = MWTSet.MWTDB.mwtpath;
%         if isequal(p1,pMWT)
%             MWTSet.timestamp = generatetimestamp;
%             MWTSetval = true;
%         end
%     end    
% end


%% create MWTSet
MWTSet.ProgramName = funame;
MWTSet.timestamp = generatetimestamp; % generate timestamp
MWTSet.PATHS.pMWT = pMWT;
MWTSet.MWTDB = parseMWTinfo(pMWT);

% record limits from varin
if numel(VarIn) > 1
    callName = VarIn(1:2:numel(VarIn));
    setValue = VarIn(2:2:numel(VarIn));
    for x = 1:numel(callName)
        MWTSet.limit.(callName{x}) = setValue{x}; 
    end
end

% make save folder
if isdir(pSave) == 0; mkdir(pSave); end
MWTSet.PATHS.pSaveA = pSave;



%% make a copy of Dance to pSave folder
if strcmp(copyfun,'on')
    ps = funpath;
    pd = sprintf('%s/%s_%s.m',pSave,funame,times);
    copyfile(ps,pd);
end











