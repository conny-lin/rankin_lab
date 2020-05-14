function D2 = import_trv_table2(pMWT,MWTDB,varargin)
%% Versions
% S = import_trv2(pMWT,varargin)

% -------------------------------------------------------------------------
%                           PROCESS INPUTS
% -------------------------------------------------------------------------
% input variables
% -------------------------------------------------------------------------
if nargin==1
    p = '/Volumes/COBOLT/MWT/MWTDB.mat';
    if ~exist(p,'file')
        error('need MWTDB input');
    else
        load(p,'MWTDB');
    end
end
legend_output = {'mwtid','tap','time','Ntrack','RevFreq', 'RevSpeed', 'RevDur'};
mwtdbvar = {'mwtid','expname','groupname','strain','rx'};
displayopt = false;
% -------------------------------------------------------------------------
%% varargin
% -------------------------------------------------------------------------
vararginProcessor
% -------------------------------------------------------------------------
%% reporting
% -------------------------------------------------------------------------
if displayopt; fprintf('starting %s\n',mfilename); end
% -------------------------------------------------------------------------
%% import into structural array
% -------------------------------------------------------------------------
load([fileparts(mfilename('fullpath')),'/legend_trv.mat'],'legend_trv');

Data = cell(size(pMWT));
for mwti =1:numel(pMWT)

    % get trv path
    pf = getpath2chorfile(pMWT(mwti),'*.trv','reporting',displayopt);
    pfval = true;
    if numel(pf)==0 % if no trv, try chor
        warning('no trv exist');
        [~,pMWTcS] = chormaster5('BeethovenOnly',pMWT(mwti));
        if ~isempty(pMWTcS)
           pf = getpath2chorfile(pMWT(mwti),'*.trv','reporting',displayopt);
        else
            pfval = false;
        end
    elseif numel(pf{1}) < numel(pMWT{mwti})
        if iscell(pf{1})
           pf1 = pf{1};
           [~,fn] = cellfun(@fileparts,pf1,'UniformOutput',0);
           fnlist = dircontent(pMWT{mwti});
           ncatch = nan(size(fn));
           for x = 1:numel(pf1)
                ncatch(x) = sum(regexpcellout(fnlist,fn{x}));
           end
           pf = pf1(ncatch==max(ncatch));
        end
    end
    

    %% see version of trv
    if pfval
        pf = char(pf);
        % read trv
        fileID = fopen(pf,'r');
        d = textscan(fileID,'%s', 2-1,'Delimiter', '', 'WhiteSpace', '');
        fclose(fileID);
        if strcmp(d{1}{1}(1),'#') ==1 % if trv file is made by Beethoven
            d = dlmread(pf,' ',5,0); 
        else % if trv file is made by Dance
            d = dlmread(pf,' ',0,0);
        end
        
        % catch error
        if size(d,2)==27
            d(:,[2 6 7 11 17 18 22]) = [];
        end
        if size(d,2) ~= numel(legend_trv); error('trv col number wrong'); end

        % add tap
        tapN = (1:size(d,1))';
        mwtid = repmat(mwti,size(d,1),1);
        time = d(:,strcmp(legend_trv,'time'));
        
        nrev = d(:,strcmp(legend_trv,'N_Rev'));
        nForPause = d(:,strcmp(legend_trv,'N_ForwardOrPause'));
        nTotal = nrev+nForPause+d(:,strcmp(legend_trv,'N_alreadyRev'));
        revfq = nrev./(nForPause + nrev);
        revdis = d(:,strcmp(legend_trv,'RevDis'));
        revdur = d(:,strcmp(legend_trv,'RevDur'));
        revspeed = revdis./revdur;
        % correct nan
        i = find(isnan(revspeed));
        revspeed(i(revdis(i)==0)) = 0;
        
        D = [mwtid tapN time nTotal revfq revspeed revdur];
        Data{mwti} = D;
        
    end
end
% -------------------------------------------------------------------------
%%  remove empty
% -------------------------------------------------------------------------
Data(cellfun(@isempty,Data),:) = [];
% -------------------------------------------------------------------------
%                   MAKE TABLE | 20170822
% -------------------------------------------------------------------------
D2 = cell2mat(Data);
D2 = array2table(D2,'VariableNames',legend_output);
%% add MWTDB variables | 20170822
D2 = outerjoin(MWTDB(:,mwtdbvar),D2,'type','right','MergeKeys',1);









