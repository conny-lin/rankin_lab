function DataMeta = IS_getData_Trinity(pMWT,mwtID,startTime,endTime)

    load('/Users/connylin/Code/proj/rankin_lab/Modules/Chor_output_handling/legend_trinity.mat')
    legend_trinity(end+1) = {'id'};
    DataMeta = cell(size(pMWT)); % create array
    mwtIDheader = cell(size(pMWT)); % create mwtid header


    for mwti = 1:numel(pMWT)
        processIntervalReporter(numel(pMWT),20,'MWT',mwti);
        pmat = char(getpath2chorfile(pMWT(mwti),'trinitySummary.mat','reporting',0)); % get pmat
        D = load(pmat); % load data

        % extract time
        D1 = D.masterData(:,2);
        nWorms = size(D1,1);
        ti = nan(nWorms,1);
        tf = ti;
        for wrmid = 1:nWorms
            ti(wrmid) = D1{wrmid,1}(1,1);
            tf(wrmid) = D1{wrmid,1}(end,1);
        end
        D = D.masterData(ti <= startTime & tf >= endTime,:); % get data within time frame

        if ~isempty(D) 

            % add worm id to last column
            ci = strcmp(legend_trinity,'id');
            a = D(:,1);
            b = regexprep(a,'^0{1,}','');
            id = cellfun(@str2num,b); % worm id
            nWorms = size(D,1);        
            for wrmi = 1:nWorms
                id1 = id(wrmi);
                n = size(D{wrmi,2},1);
                D{wrmi,2}(:,ci) = repmat(id1,n,1);
            end 

            % convert other datadata
            D = cell2mat(D(:,2)); 
            t = D(:,strcmp(legend_trinity,'time'));
            D(~(t >= startTime & t <= endTime),:) = [];

            % calculate speed, speedbm, curve by worms
            speed = D(:,strcmp(legend_trinity,'speed'));
            id = D(:,strcmp(legend_trinity,'id'));
            midline = D(:,strcmp(legend_trinity,'midline'));
            curve = D(:,strcmp(legend_trinity,'curve'));
            speedbm = speed./midline;

            clear S
            S = grpstatsTable(speed,id,'gnameutitle','wrmid');
            S.wrmid = cellfun(@str2num,S.wrmid);
            SB = grpstatsTable(speedbm,id,'gnameutitle','wrmid');
            SB.wrmid = cellfun(@str2num,SB.wrmid);
            C = grpstatsTable(curve,id,'gnameutitle','wrmid');
            C.wrmid = cellfun(@str2num,C.wrmid);

            if ~isequal(S.wrmid, SB.wrmid) || ~isequal(SB.wrmid,C.wrmid)
                error('wormid does not match');
            end
            D1 = [S.wrmid C.mean SB.mean S.mean];

            % put in main datasheet
            DataMeta{mwti} = D1;
            mwtIDheader{mwti} = repmat(mwtID(mwti),size(D1,1),1);

        end
    end
    % collapse data
    Data = [cell2mat(mwtIDheader) cell2mat(DataMeta)];
    DataMeta = array2table(Data,'VariableNames',{'mwtid','wrmid','curve','speedbm','speed'});
end