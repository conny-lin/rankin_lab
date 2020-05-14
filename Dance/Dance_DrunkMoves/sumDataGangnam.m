function [Data,MWTDB] = sumDataGangnam(pMWT,timeset,varargin)

% settings: +++++++
VarName = {'speed','curve'};
vararginProcessor
% ------------------

% convert settings +++++++
VNAMES = [{'id','time'} VarName];

tstart = timeset(1);
tend = timeset(3);
tint = timeset(2);
timestartlist = [tstart:tint:tend];
% ----------------------

%% Chor Gangnam ============================================================
% find existing gangnam
[~,~,pMWTmiss] = getpath2chorfile(pMWT,'Gangnam.mat');
% chor plates without gangnam
[~,~,pMWTS] = getpath2chorfile(pMWTmiss,'*.gangnam.*.dat');
[Legend,~,~,~] = chormaster4('Gangnam',pMWTS);
legend_gangnam = Legend{1}{2};
inputname = Legend{1}{1};
% convert gangnam into .mat
convertchorNall2mat(pMWT,inputname,'Gangnam');

% check if mwt plates have gangnam 
pMat = getpath2chorfile(pMWT,'Gangnam.mat');
% reconstruct MWT database
pMWT = cellfun(@fileparts,pMat,'UniformOutput',0);
MWTDB = parseMWTinfo(pMWT);
% -----------------------------------------------------
%% =========================================================================


% declare input arrays +++++++++
nMWT = numel(pMat);
sampleSize = nan(nMWT,1);
DataAll = cell(nMWT,1);
% -----------------------------


for mwti = 1:nMWT
    
    % report progress ++++++++++++++++++++
    processIntervalReporter(nMWT,10,'MWT',mwti);
    % -----------------------------------
    
    % get cycle specific info +++++++++
    pmat = pMat{mwti};
    pmwt = fileparts(pmat);
    % ---------------------------------
    
    % load data ++++++++++++++++
    D = load(pmat);
    i = D.time(:,1) <= tstart & D.time(:,2) >= tend;
    files = D.Data(i);
    % ---------------------------
    
    % get sample size +++++++++++++
    n = numel(files);
    sampleSize(mwti) = n;
    % ----------------------------
    
    if n>0
        % extract data +++++++++++++
        T_plate = table;
        for wrmi = 1:numel(files)
            % get only needed data ++++++++
            d = files{wrmi};
            [a,legend] = extractChorData(d,legend_gangnam,'varnames',VNAMES,'outputtype','array');
            a = array2table(a,'VariableNames',legend);
            % ----------------------------
            
            
            % condense into time list +++++++++
            a.timer = round(a.time);
            tu = unique(a.timer);
            tun = numel(tu);
            T = table;
            T.time = tu;
            T.id = repmat(unique(a.id),tun,1);
            
            VarName = {'speed','curve'};
            for si = 1:numel(VarName)
                legname = VarName{si};
                T2 =  statsBasicG(a.(legname), a.timer);
                
                T3 = table;
                T3.time = T2.gname;
                T3.(legname) = T2.mean;
                T = outerjoin(T,T3,'MergeKeys',1,'Type','left');                
            end
            % -----------------------------------
            
            % add to plate table +++++
            T_plate = [T_plate; T];
            % -----------------------
            
        end
        
        
        % calculate plate stats
        
        T = table;
        T.time = timestartlist';
        
        for si = 1:numel(VarName)
            legname = VarName{si};
            T2 =  statsBasicG(T_plate.(legname), T_plate.time);
            
            T3 = table;
            T3.time = T2.gname;
            T3.(legname) = T2.mean;
            
            T = outerjoin(T,T3,'MergeKeys',1,'Type','left'); 
        end
        
        
        % create mwt paths for each row ++++++++
        leg = table;
        leg.mwtpath = repmat(mwti, size(T,1),1);
        % ---------------------------------------
        
        % add to central table +++++++++
        T = [leg T];
        T = table2array(T);
        DataAll{mwti} = T;
        % ---------------------------
        
    end
    
end

% convert to table +++++++
DataAll(cellfun(@isempty,DataAll)) = [];
a = cell2mat(DataAll);
Data = array2table(a,'VariableNames',[{'mwtid','time'}, VarName]);