function [D,Legend] = import_drunkposture2(pMWT,MDb,Vind)
%% import_drunkposture2(pMWT) 


%% IMPORT RAW DATA & COMPILE TO .MAT (.drunkposture2.dat)
Legend = {'time' 'number' 'goodnumber' 'speed' 'length' 'width' 'aspect',...
    'kink' 'bias' 'curve' 'area' 'midline' 'morphwidth'}';
% get info needed
importname = 'drunkposture2';

VNames = fieldnames(Vind);

% code (only import the first file name)
fprintf('Importing %s...\n',importname); 
D = table;
val = false(size(pMWT));
for mwti = 1:numel(pMWT);
    [~,fn] = dircontent(pMWT{mwti},'*.drunkposture2.dat');  
    % delete fntemp
    fntemp = fn(regexpcellout(fn,'\<[.]'));
    fn(regexpcellout(fn,'\<[.]')) = [];
    for x = 1:numel(fntemp)
        fprintf('- delete temp files %s\n',fntemp{x})
        delete(fntemp{x});
    end
    if isempty(fn) == 0 % get data, if there is a chor output
        val(mwti) = true;
        fn = fn{1};
        fnsave = regexprep(fn,'[.]dat\>','.mat');
        [~,fnmat] = fileparts(fnsave);
        [~,fn2] = dircontent(pMWT{mwti},[fnmat,'.mat']);
        if isempty(fn2) == 0
            d = load(fn2{1});
            d = d.d;
        else
            d = dlmread(fn);
            d = array2table(d,'VariableNames',Legend);
            save(fnsave,'d');
        end
        % get rows of d
        nRow = size(d,1);
        [~,fmwt] = fileparts(pMWT{mwti});
        % get database index
        info = MDb(ismember(MDb.mwtname,fmwt),:);
        % replace info with varindex
        for x = 1:numel(VNames)
            ind = find(ismember(Vind.(VNames{x}),info.(VNames{x})));
            info.(VNames{x}) = ind;
        end
        info = repmat(info,nRow,1);
        recordseqN = array2table((1:nRow)','VariableNames',{'recordseq'});
        d = [info,recordseqN,d];
        % combine data
        D = [D;d];
    end
end
% clean up and report
% remove files without import
% M = MDb;
% M.(['import_',importname]) = val;
% MWTSet.Info.MWTDb_status = M;
% MWTSet.Info.MWTDbInd_original = MWTSet.Info.MWTDbInd;
% MWTSet.Info.MWTDbInd(~val,:) = [];
% MWTSet.Info.MWTDb_original = MWTSet.Info.MWTDb;
% MWTSet.Info.MWTDb(~val,:) = [];
% fprintf('-%d/%d MWT does not have chor output: removed\n',sum(val == false),numel(i));

% store results
% MWTSet.Import.(importname) = D;
fprintf('Finish importing [%s]\n',importname);
