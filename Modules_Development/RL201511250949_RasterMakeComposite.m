%% for 60sISI

%% Paths
pAnalysis = '/Volumes/COBOLT/MWT';
pData = '/Volumes/COBOLT/MWT';
% add function
addpath('/Users/connylin/Dropbox/MATLAB/Function_Library_Public');

%% input variables 60sISI
preplateTime = 100;
tStart = 85;
isiT = 10;
tEnd = (isiT*30)+100+10;
tInt = 10;
timePairs = [tStart tStart+tInt:isiT:tEnd-tInt; tStart+tInt tStart+tInt*2:isiT:tEnd]';
tapNumber = 30;

%% get database
addpath('/Users/connylin/Dropbox/RL/Code/Library/Modules/MWTDatabase');
MWTD = load(sprintf('%s/MWTDatabase.mat',pData));
MWTD = MWTD.MWTDatabase;

% find liquid transfer date
a = find(regexpcellout(MWTD.mwt.groupname,'Liquid'));
liquidTsfDate = MWTD.mwt.exp_date(a(end));

expsetNameList = MWTD.expset.expset;
for expsetNameListi = 11%:numel(expsetNameList)
    expsetName = expsetNameList{expsetNameListi};
    
    % create experiment set folder
    pDest = sprintf('/Users/connylin/Dropbox/RL/MWT_Analysis_ByExp/%s %dISI',expsetName,isiT);
    if isdir(pDest) == 0; mkdir(pDest); end
    
    % find genes from exp set
    genes = MWTD.expset.genes{ismember(MWTD.expset.expset,expsetName)};
    val = false(size(MWTD.mwt,1),1);
    for x = 1:numel(genes)
        val(~cellfun(@isempty,regexp(MWTD.mwt.genotype,genes(x)))) = true;
    end
    strainT = MWTD.mwt.strain(val);
    
    % find targets
    A = MWTD.mwt;
    i = ismember(A.strain,strainT) &...
        ismember(A.rx,{'NA','400mM'}) &...
        A.ISI == isiT &...
        A.preplate == preplateTime &...
        A.tapN == tapNumber &...
        A.exp_date > liquidTsfDate;
    if sum(i) > 0; 
        D = MWTD.mwt(i,:);

        % find experiments with controls
        expT = unique(D.expname);
        i = ismember(A.expname,expT) & ...
            ismember(A.strain,{'N2'}) & ...
            ismember(A.rx,{'NA','400mM'}) &...
            A.ISI == isiT &...
            A.preplate == preplateTime &...
            A.tapN == tapNumber &...
            A.exp_date > liquidTsfDate;
        D = [D;A(i,:)];

        % eliminate duplicates
        [~,i] = unique(D.mwt_id);
        D(~i,:) = [];

        % display names
        display('group names:');
        disp(char(unique(D.groupname)));

        % create outputs
        writetable(D,[pDest,'/plate_info.csv']);
        % make table
        save([pDest,'/plate_info.mat'],'D');


        %% chor Trinity
        % look if all files have required files
        reqFile = 'trinitySummary.mat';
        pMWT = D.mwtpath;
        val = false(size(pMWT));
        for mwti = 1:numel(pMWT)
            a = dircontent(pMWT{mwti},reqFile);
            if isempty(a) == 0; val(mwti) = true; end
        end
        fprintf('%d/%d MWT do not have %s file\n',sum(~val), numel(val),reqFile);

        % check if no reqFile has summary file
        reqFile = '*trinity*.dat';
        pMWT = pMWT(~val);
        val = false(size(pMWT));
        for mwti = 1:numel(pMWT)
            a = dircontent(pMWT{mwti},reqFile);
            if isempty(a) == 0; val(mwti) = true; end
        end
        fprintf('%d/%d MWT have %s file\n',sum(val), numel(val),reqFile);

        % chor
        pMWT = pMWT(~val);
        addpath('/Users/connylin/Dropbox/RL/Code/Library/Modules/Chor');
        chormaster4('Trinity',pMWT);



        %% run raster plot (use pMWT)
        % look if all files have required files
        reqFile = 'trinitySummary.mat';
        pMWT = D.mwtpath;
        val = false(size(pMWT));
        for mwti = 1:numel(pMWT)
            a = dircontent(pMWT{mwti},'trinitySummary.mat');
            if isempty(a) == 0; 
                val(mwti) = true; 
            else
                a = dircontent(pMWT{mwti},'*trinity*.dat');
                if isempty(a) == 0
                    val(mwti) = true; 
                end
            end
        end
        fprintf('%d/%d MWT do not have %s file\n',sum(~val), numel(val),reqFile);
        D(~val,:) = [];
        % create outputs
        writetable(D,[pDest,'/plate_info.csv']);
        % make table
        save([pDest,'/plate_info.mat'],'D');
        % group paths into groups
        fG = D.groupname;
        fGU = unique(D.groupname);
        pathGrouped = cell(size(fGU));
        for x = 1:numel(fGU)
            pathGrouped{x} = D.mwtpath(ismember(fG,fGU{x}));
        end

        % run
        addpath('/Users/connylin/Dropbox/RL/Code/Library/Modules/rasterPlot_colorSpeed');

        for gi = 1:numel(fGU)
            fprintf('Raster plotting group [%s]...\n',fGU{gi});
            for ti = 1:size(timePairs,1)
                % reporting
                str = sprintf('rasterPlot_%d_%d*',timePairs(ti,1),timePairs(ti,2));
                fprintf('--raster plot [%d/%d]\n',ti,size(timePairs,1));
                %% run raster plot
                [~,~,~,Tog] = rasterPlot_colorSpeed(pathGrouped{gi},timePairs(ti,1),timePairs(ti,2),'NWORMS',Inf,'InputType',2,'visibleG',0);
                TOG{ti} = Tog;
                
                %%
                a1 = subplot(2,1,1)
                imagesc(Tog,'Parent',a1)
                
                
                
                
                
                return
                % create save folder
                pSave = [pDest,'/',fGU{gi},'/rasterPlot_colorSpeed'];
                if isdir(pSave) == 0; mkdir(pSave); end
                cd(pSave);
                set(figure1,'PaperPositionMode','auto'); % set to save as appeared on screen
                print (figure1,'-depsc', '-r1200', savename); % save as eps
                close;
            end
        end


        % report done
        fprintf('done[ %s]*\n',expsetName);
    end
end