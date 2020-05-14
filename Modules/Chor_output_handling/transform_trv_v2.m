function [D,PF] = transform_trv_v2(Data,varargin)
%% transform_trv_v2(A)
% A is a cell array containing trv imports
% Data = MWTSet.Import.trv;
 

%% DEFAULTS AND VARARGIN
tapindx = 1; % will create tap index
% process varargin
vararginProcessor; 

%% MAKING SENSE OF TRV (r20151126)
% get data
A = Data.data;
% get path
pMWT = Data.mwtpath;

%% parse info from pMWT
% Get these info:
% fn = {'mwtpath',GroupVarName};
% Prefix = table;
% for x = 1:numel(fn)
%    Prefix.(fn{x}) = MWTDBInd.(fn{x});
% end

%% get essential variables
var = {'time' 'N_alreadyRev' 'N_ForwardOrPause' 'N_Rev' 'RevDis' 'RevDur' 'RevDur_SE'};
D  = table;
PF = {};
for mwti = 1:size(Data,1)
    % get data
    d = A{mwti};
    if isempty(d) == 1
       error('no data found'); 
    else % get only essential var
        a = table;
        for vi = 1:numel(var); a.(var{vi}) = d.(var{vi}); end
        % make tap time index
        if tapindx == 1
           t = table;
           t.tap = (1:size(a,1))';
           a = [t a];
        end
        
    end
    % make prefix
    pfx = repmat(pMWT(mwti,:), size(d,1),1);
    % add prefix to data   
%     d = [pfx a];
    PF = [PF;pfx];
    % add to master table;
    D = [D;a];
end


%% calculate reversal frequency 
% frequency is # reversed divided by total number not already reversing
D.RevFreq = D.N_Rev./(D.N_Rev + D.N_ForwardOrPause);

% calculate reversal speed
% reversal speed is reversal distance divided by reversal duration
D.RevSpeed = D.RevDis./D.RevDur;




