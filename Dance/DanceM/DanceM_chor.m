function [MWTSet,MWTDB] = DanceM_chor(MWTSet,chortype,varargin)

%% default
displayopt = false;
overwrite = false;
%% varargin
vararginProcessor;

pMWT = MWTSet.pMWT;

%% run chor
[Legend,pMWTpass,pMWTfailed] = chormaster5(chortype,pMWT,'deletedat',1,'displayopt',displayopt); % chor

%% post chor modification
[MWTSet,~] = DanceM_MWTSet_postChor(MWTSet,pMWT,pMWTpass,pMWTfailed,Legend);
MWTDB = MWTSet.MWTDB;

%% patches for post data transformation
choroutputnames = fieldnames(Legend);
for ci = 1:numel(choroutputnames)
    cn = choroutputnames{ci};
        switch cn
            case 'trinity'
                for i = 1:numel(pMWTpass)
                    if displayopt
                        loopreporter(i,'mwt',1,numel(pMWT)); 
                    end
                    revert_trinitySummary2indivialFiles(pMWTpass{i},overwrite,'displayopt',displayopt)
                end
            case 'drunkposture2'
                MWTDB = parseToleranceName(MWTDB);
                MWTSet.MWTDB = MWTDB;
        end
end


