function [pFun] = DanceM_funpath(packagepath)
%% DanceM_funpath
% add paths to general functions & define general paths
% e.g. packagepath = '/Users/connylin/Dropbox/rl/Code/Modules/Dance/Dance_DrunkMoves';


%% create common paths
pFun = struct;
pFun.pMWTData = '/Volumes/COBOLT/MWT';
pFun.pMWTDB = '/Users/connylin/Dropbox/RL/MWTDB';


%% matlab library
pMatFun = '/Users/connylin/Dropbox/Code/Matlab/Library';
% add packges
addpath([pMatFun,'/General']);
addpath([pMatFun,'/Stats']);
addpath([pMatFun,'/Graphs']);

%% rankin lab library
pRLFun = '/Users/connylin/Dropbox/RL/Code/Modules';
% add packges
addpath([pRLFun,'/Chor']);
addpath([pRLFun,'/Graphs']);
addpath([pRLFun,'/Graphs/rasterPlot_colorSpeed']);
addpath([pRLFun,'/MWTDatabase']);

%% add Dance package modules if any
if nargin > 0
    p = fileparts(packagepath);
    if isdir(p) == 0; error('path given is not a folder'); end
    [~,~,~,pf] = dircontent(p);
    for x = 1:numel(pf)
        addpath(pf{x});
    end
end

