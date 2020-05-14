function [pMWTp] = convertTrinityDat2Mat(pMWT,deleteDat)


if nargin < 2 
    deleteDat = 1;
end

%% transform data into mat file
trinitymissing = false(size(pMWT));
for z = 1:numel(pMWT)
    fprintf('%d/%d pMWT\n',z,numel(pMWT)); % progress report
    [~,TrinityD] = dircontent(pMWT{z},'trinitySummary.mat');
    
    if numel(TrinityD) == 0
        [~,pTrinity] = dircontent(pMWT{z},'*.trinity.*.dat');
        i = regexpcellout(pTrinity,'\<[.]*');
        cellfun(@delete,pTrinity(i));
        pTrinity(i) = [];
        
        if isempty(pTrinity) == 0
            masterData = cell(size(pTrinity,1),2); % declare output array
            [~,f] = cellfun(@fileparts,pTrinity,'UniformOutput',0); % get worm ID
            masterData(:,1) = regexpcellout(f,'(?<=trinity[.])\d{1,}','match'); % record worm id
            for wormi = 1:numel(pTrinity)
                masterData{wormi,2} = dlmread(pTrinity{wormi}); % read .trinity.*.dat data from a worm
            end
            save([pMWT{z},'/trinitySummary.mat'],'masterData','-v7.3'); % save file as .mat
        else
            trinitymissing(z) = true;
        end
    end
    
    % clear trinity.dat
    [~,TrinityD] = dircontent(pMWT{z},'trinitySummary.mat');
    if numel(TrinityD) ==1
        [~,pTrinity] = dircontent(pMWT{z},'*.trinity.*.dat');
        if deleteDat && numel(pTrinity) > 1
            cellfun(@delete,pTrinity); 
        end 
    end

end

pMWTp = pMWT(trinitymissing);


