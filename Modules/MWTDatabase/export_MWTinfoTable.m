function T = export_MWTinfoTable(pMWT,varargin)

%% DEFAULTS AND PROCESS VARARGIN
pSave = '/Users/connylin/Dropbox/RL/Dance Output'; % default save output path
suffix = ''; % suffix of save file name
display = false; % default to display table result
savetable = true; % defaut to not display
% process varargin
vararginProcessor;

t = parseMWTinfo(pMWT);
T = table;
T.mwtpath = t.mwtpath;
T.expname = t.expname;
T.groupname = t.groupname;
T.mwtname = t.mwtname;

%% OPTIONAL
% add suffix (optional)
if isempty(suffix) == 1
    savename = sprintf('%s/mwt_summary.csv',pSave);
else
    savename = sprintf('%s/mwt_summary %s.csv',pSave,suffix);
end

% save (optional)
if savetable == true
    writetable(T,savename);
end

% display (optional)
if display == true; disp(t); end


 
end