function Trinity = load_Trinity_v1707(pMWTp)

%% find and load trinity summary file, optional clean up
Trinity = [];

pTrinityMat = dircontent(pMWTp,'trinitySummary.mat');

if isempty(pTrinityMat) % if no mat files
   if displayopt; fprintf('no trinitySummary.mat file found\n'); end
   [~,p] = dircontent(pMWTp,'*trinity.*.dat'); % look for dat files
   % check if trinity.dat exist
   if isempty(p)  % if no .dat
       return
   else % if has .dat, convert
       convertTrinityDat2Mat({pMWTp},1); % generate trinity.mat
   end        
end
% load trinity file
pTrinityMat = fullfile(pMWTp,'trinitySummary.mat');
Trinity = load(pTrinityMat,'masterData'); 
Trinity = Trinity.masterData;
% -----------------------------------------------------------------


