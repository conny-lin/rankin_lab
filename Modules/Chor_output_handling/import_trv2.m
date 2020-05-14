function S = import_trv2(pMWT,varargin)
%% Versions
% S = import_trv2(pMWT,varargin)


%% input variables
% get necessary legends

legend_output = {'tap','time', 'N_alreadyRev', 'N_ForwardOrPause', 'N_Rev', 'RevDis', 'RevDur'};
vararginProcessor


%% reporting
fprintf('starting %s\n',mfilename);


%% load variables
pTrvLegend = [fileparts(mfilename('fullpath')),'/legend_trv.mat'];
load(pTrvLegend);
% process input variables
ind_get = ismember(legend_trv,legend_output);


%% import into structural array
S = struct;
for mwti =1:numel(pMWT)
    % prepare summary output identifiers
    db = parseMWTinfo(pMWT{mwti});
    S(mwti).mwtpath = pMWT{mwti};
    S(mwti).expname = db.expname;
    S(mwti).groupname = db.groupname;
    S(mwti).mwtname =  db.mwtname;
    % get trv path
    pf = getpath2chorfile(pMWT(mwti),'*.trv','reporting',0);
    if numel(pf) > 1; 
        error('more than one trv paths'); 
    elseif numel(pf)==0
        warning('no trv exist');
        [~,pMWTcS] = chormaster5('BeethovenOnly',pMWT(mwti));
        if ~isempty(pMWTcS)
           pf = getpath2chorfile(pMWT(mwti),'*.trv','reporting',0);
        else
            pfval = false;
        end
    else
        pfval = true;
    end
    
    % see version of trv
    if pfval
        d = dlmread(char(pf));
        ncol = size(d,2);
        if ncol ~= numel(legend_trv); error('trv col number wrong'); end
        % add tap
        tapnumber = (1:size(d,1))';
        D = array2table([tapnumber,d(:,ind_get)],'VariableNames',legend_output);
        % frequency is # reversed divided by total number not already reversing
        D.RevFreq = D.N_Rev./(D.N_Rev + D.N_ForwardOrPause);
        % calculate reversal speed: reversal speed is reversal distance divided by reversal duration
        D.RevSpeed = D.RevDis./D.RevDur;
        D.RevDis = [];

        % enter variable in structural array
        vnames = D.Properties.VariableNames;
        for vi = 1:numel(vnames)
            vn = vnames{vi};
            S(mwti).(vn) = D.(vn);
        end
    end
end






