function [GroupedpMWT4Dance] = groupMWTpaths4Dance(pMWT,StrainNameTarget)

%% parse pMWT into 
T = parseMWTinfo(pMWT);


%% check if there is a strain name target
if ~isvarname('StrainNameTarget')    
    StrainNameTarget = T.strain;
    % get rid of wildtype
    StrainNameTarget(ismember(StrainNameTarget,'N2')) = [];
end

%% group pMWT with N2, N2_400mM, strain, strain_400mM

s = T.strain;
GroupedpMWT4Dance = cell(numel(StrainNameTarget),1);
for a = 1:numel(StrainNameTarget)
    
    b = T.mwtpath(ismember(s,{'N2',StrainNameTarget{a}}));
    GroupedpMWT4Dance{a} = b;
     
end


end