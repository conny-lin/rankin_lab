function DataMeta = IS_getData(pMWT,mwtID,startTime,endTime)

    load('/Users/connylin/Dropbox/Code/Matlab/Library RL/Modules/Chor_output_handling/legend_gangnam.mat')

DataMeta = cell(size(pMWT));
mwtIDheader = cell(size(pMWT));

for mwti = 1:numel(pMWT)
    processIntervalReporter(numel(pMWT),20,'MWT',mwti);
    pmat = char(getpath2chorfile(pMWT(mwti),'Gangnam.mat','reporting',0));
    D = load(pmat);
    D = D.Data(D.time(:,1) <= startTime & D.time(:,2) >= endTime);
    if ~isempty(D); 
        D = cell2mat(D); 
        t = D(:,strcmp(legend_gangnam,'time'));
        D(~(t >= startTime & t <= endTime),:) = [];
        
        % calculate speed, speedbm, curve by worms
        speed = D(:,strcmp(legend_gangnam,'speed'));
        id = D(:,strcmp(legend_gangnam,'id'));
        midline = D(:,strcmp(legend_gangnam,'midline'));
        curve = D(:,strcmp(legend_gangnam,'curve'));
        speedbm = speed./midline;

        S = grpstatsTable(speed,id,'gnameutitle','wrmid');
        SB = grpstatsTable(speedbm,id,'gnameutitle','wrmid');
        C = grpstatsTable(curve,id,'gnameutitle','wrmid');
        if ~isequal(S.wrmid, SB.wrmid) || ~isequal(SB.wrmid,C.wrmid)
            error('wormid does not match');
        end
        D = [S.wrmid C.mean SB.mean S.mean];

        % put in main datasheet
        DataMeta{mwti} = D;
        mwtIDheader{mwti} = repmat(mwtID(mwti),size(D,1),1);
    end
end
% collapse data
Data = [cell2mat(mwtIDheader) cell2mat(DataMeta)];
DataMeta = array2table(Data,'VariableNames',{'mwtid','wrmid','curve','speedbm','speed'});

