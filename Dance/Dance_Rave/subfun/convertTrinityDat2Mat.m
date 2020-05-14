function convertTrinityDat2Mat(pMWT,InputType)

%% transform data into mat file
for z = 1:numel(pMWT)
    [~,TrinityD] = dircontent(pMWT{z},'trinitySummary.mat');
    if numel(TrinityD) == 0
        [~,pTrinity] = dircontent(pMWT{z},'*.trinity.*.dat');
        i = regexpcellout(pTrinity,'\<[.]*');
        cellfun(@delete,pTrinity(i));
        pTrinity(i) = [];
        if isempty(pTrinity) == 0
            masterData = cell(size(pTrinity,1),2); % declare output array
            % get worm ID
            [~,f] = cellfun(@fileparts,pTrinity,'UniformOutput',0);
            % record worm id
            masterData(:,1) = regexpcellout(f,'(?<=trinity[.])\d{1,}','match');
            for wormi = 1:numel(pTrinity);
                masterData{wormi,2} = dlmread(pTrinity{wormi}); % read .trinity.*.dat data from a worm
            end
            save([pMWT{z},'/trinitySummary.mat'],'masterData','-v7.3'); % save file as .mat

            if InputType ==1
                cellfun(@delete,pTrinity); % clear trinity.dat
            end
        end
    end
end