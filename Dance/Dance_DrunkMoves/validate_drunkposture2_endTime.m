function [D,Dbad] = validate_drunkposture2_endTime(D,pMWT,Vind)
% end time has to be the same as what's said on experiment
% D = Import;
% calculate expected tap time
A = parseMWTinfo(pMWT);
t = A.preplate + (A.ISI.*A.tapN) + A.postrec;
% calculate max time per plate
[gn,mx] = grpstats(D.time,D.mwtname,{'gname','max'});
% see which plate has bad end time
mwtname = Vind.mwtname(cellfun(@str2num,gn));
[i,j] = ismember(mwtname,A.mwtname);
i = (t(j(i)) - mx) > 0;
badmwt = mwtname(i);
[i,j] = ismember(badmwt,Vind.mwtname);
badmwtid = j(i);
i = ismember(D.mwtname,badmwtid);
Dbad = D(i,:);
D(i,:) = [];