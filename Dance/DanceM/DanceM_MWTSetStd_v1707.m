function [MWTSet,pSaveA] = DanceM_MWTSetStd_v1707(pMWT,funpath,VarIn,pSave)
% DanceM_MWTSetStd2 - use DanceM_createIndex2 for flexibility

%% default
% useoldfile = false;
% timestamp = 'off';
% outputsuffix = '';
% copyfun = 'on';
% vararginProcessor;


%% save folder
funpath = [funpath,'.m'];
[~,funame] = fileparts(funpath);
pSaveA = create_savefolder(pSave,funame);
times = generatetimestamp;


%% create MWTSet
MWTSet.ProgramName = funame;
MWTSet.timestamp = times; % generate timestamp
MWTSet.pSave = pSaveA; % save folder
MWTSet.pMWT = pMWT;
MWTSet.MWTDB = parseMWTinfo(pMWT);


%% record limits from varin
if numel(VarIn) > 1
    callName = VarIn(1:2:numel(VarIn));
    setValue = VarIn(2:2:numel(VarIn));
    for x = 1:numel(callName)
        MWTSet.setting.(callName{x}) = setValue{x}; 
    end
end



%% make a copy of Dance to pSave folder
% if strcmp(copyfun,'on')
ps = funpath;
pds = create_savefolder(pSaveA,'_code');
pd = fullfile(pds,sprintf('%s_%s.m',funame,times));
copyfile(ps,pd);
% end











