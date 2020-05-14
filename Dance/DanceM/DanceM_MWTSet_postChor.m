function [MWTSet,pMWTpass] = DanceM_MWTSet_postChor(MWTSet,pMWT,pMWTpass,pMWTfailed,Legend)

% update MWTDB 
MWTSet.MWTDB_input = parseMWTinfo(pMWT);
MWTSet.MWTDB_failed = parseMWTinfo(pMWTfailed);
MWTSet.MWTDB = parseMWTinfo(pMWTpass);
MWTSet.legend_chor = Legend;