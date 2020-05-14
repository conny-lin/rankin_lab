% MWT output info sheet
a = MWTSet.Info.MWTDB;
writetable(a,sprintf('%s/info.csv',pSave));
% save MWTSet
cd(pSave);
save(mfilename,'MWTSet');

%% FINISH
fprintf('*** completed ***\n');