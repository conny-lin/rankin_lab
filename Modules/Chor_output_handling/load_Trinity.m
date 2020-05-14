function Trinity = load_Trinity(pMWTp)

cleanup = 0;
%% find and load trinity summary file, optional clean up
Trinity = [];
if isempty(dircontent(pMWTp,'trinitySummary.mat'))
%     D = load([pMWTp,'/trinitySummary.mat'],'masterData'); 
%     D = D.masterData;
% else
%    fprintf('no trinitySummary.mat file found\n');
   [~,p] = dircontent(pMWTp,'*trinity.*.dat'); 
   if isempty(p) == 0
        convertTrinityDat2Mat({pMWTp});
   else
       warning('no trinity.*.dat file found, need to chor, skip');
       return
   end
end
D = load([pMWTp,'/trinitySummary.mat'],'masterData'); 
D = D.masterData;
Trinity = D;


if cleanup == 1
   [~,p] = dircontent(pMWTp,'*trinity.*.dat'); 
   if isempty(p) == 0
       fprintf('cleaning up trinity.*.dat file\n');
       cellfun(@delete,p);
   end
end

