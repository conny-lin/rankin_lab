function pMWTval = chor_Trinity(pMWT)

%%


pMWTc = convertTrinityDat2Mat(pMWT,1); 
L = chormaster4('Trinity',pMWTc);
% summarize trinity data and delete .dat file to save memory
pMWTbad = convertTrinityDat2Mat(pMWTc,1); 
% exclude bad files
pMWToriginal = pMWT;
pMWTval = pMWT(~ismember(pMWT,pMWTbad));
