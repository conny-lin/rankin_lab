%% new raster
% - sort by speed
% - color adjusted to jet

%% initializing
clc; clear; close all;

%% paths
% function paths
addpath('/Users/connylin/Dropbox/Code/Matlab/Library/General');
addpath('/Users/connylin/Dropbox/RL/Code/Modules/Graphs/rasterPlot_colorSpeed');
% save and data paths
pData = fileparts(mfilename('fullpath'));

%% get group folders
[~,~,groupfolder, pG] = dircontent(pData);

%% loop for group folders
for gi = 1:numel(pG)
    pg = pG{gi};
    % get filenames
    filenames = dircontent(pg,'*rasterData.mat');
    a = regexpcellout(filenames,'_','split');
    a = a(:,2:3);
    savesuf = cell(size(a,1),1);
    for x = 1:size(a,1)
       savesuf(x) = {strjoin(a(x,:),'-')};
    end
    % run plots for each file names
    for x = 1:numel(filenames)
        % load data
        cd(pg); load(filenames{x});
        savesuffix = savesuf{x}; % set suffix
        time = rTime(1,:);

        % plot gradient raster plot
        figure1 = rasterPlot_colorSpeed_gradient(Data,time,0);
        pS = sprintf('%s/raster gradient',pg);
        if isdir(pS) == 0; mkdir(pS); end
        titlename = sprintf('%s/raster %s',pS,savesuffix);
        print(figure1,'-depsc', titlename); % save as (r600 will give better resolution)
        print(figure1,'-dtiff', titlename); % save as (r600 will give better resolution)

        close;

        % plot sort proportion
        figure1 = moveproportion(Data,time,0);
        pS = sprintf('%s/speed dir proportion plot',pg);
        if isdir(pS) == 0; mkdir(pS); end
        titlename = sprintf('%s/speed_dir_proportion_%s',pS,savesuffix);
        print(figure1,'-dtiff', titlename); % save as (r600 will give better resolution)
        print(figure1,'-depsc', titlename); % save as (r600 will give better resolution)

        close;

    end
end



%% assemble pictures -----------------------------------------------------
timenames = {'98-128','238-268','368-398'};
targetfigureName = {'speed dir proportion plot','raster gradient'};

[~,gname] = cellfun(@fileparts,pG,'UniformOutput',0);

% find strain sets
a = regexpcellout(gname,'_','split');
i = cellfun(@isempty,a(:,2));
strain = a(i,1);
strain_condition = gname(~i);
group_condition = [strain strain_condition];
condition = a(~i,2);
pGCI = [pG(i) pG(~i)];
for tfni = 1:numel(targetfigureName);
    figureprefix = targetfigureName{tfni};
    for gci = 1:size(group_condition,1)
        gname = group_condition(gci,:);
        IM = cell(numel(timenames),numel(gname));
        pgc = pGCI(gci,:);
        for gi = 1:numel(gname)
            pg = pgc{gi};
            [~,gn] = fileparts(pg);

            pfolder = sprintf('%s/%s',pg,figureprefix);
            [figurename,pfigure] = dircontent(pfolder);
            pfigure(~regexpcellout(pfigure,'[.]tif\>')) = [];
            % get images
            for ti = 1:numel(timenames)
                p = pfigure{regexpcellout(pfigure,timenames{ti})};
                [~,fn] = fileparts(p);
                fprintf('loading file: %s/%s\n',gn,fn);
                I = imread(p);
                % trim figure
                I = I(50:end-50,100:end-100,:);
                IM{ti,gi} = I;
            end
        end

        IMM = cell2mat(IM);

        % add header
        [nrow,ncol,nlayer] = size(IMM);
        Header = repmat(255,300,ncol,3);
        IMMH = [Header;IMM];
        [nrow,ncol,nlayer] = size(IMMH);
        SideBar = repmat(255,nrow,200,3);
        IMMH = [SideBar IMMH];
        figure1 = figure('Visible','off');
        imshow(IMMH)
        % y axis
        fs = 15; leftalign = 185;
        text(leftalign,700,'% total','Rotation',90,'FontSize',fs,'HorizontalAlignment','center','VerticalAlignment','middle');
        text(leftalign,1500,'% total','Rotation',90,'FontSize',fs,'HorizontalAlignment','center','VerticalAlignment','middle');
        text(leftalign,2300,'% total','Rotation',90,'FontSize',fs,'HorizontalAlignment','center','VerticalAlignment','middle');
        % time names
        fs = 20;     leftalign = 80;     
        text(leftalign,700,timenames{1},'Rotation',90,'FontSize',fs,'HorizontalAlignment','center','VerticalAlignment','middle');
        text(leftalign,1500,timenames{2},'Rotation',90,'FontSize',fs,'HorizontalAlignment','center','VerticalAlignment','middle');
        text(leftalign,2300,timenames{3},'Rotation',90,'FontSize',fs,'HorizontalAlignment','center','VerticalAlignment','middle');

        % concentration
        topalign = 220; fs = 20;
        text(750,topalign,'0mM','Rotation',0,'FontSize',fs,'HorizontalAlignment','center','VerticalAlignment','middle');
        text(1700,topalign,condition{gci},'Rotation',0,'FontSize',fs,'HorizontalAlignment','center','VerticalAlignment','middle');
        % strain
        topalign = 100; fs = 20;
        text(1200,topalign,gname{1},'Rotation',0,'FontSize',fs,'HorizontalAlignment','center','VerticalAlignment','middle');

        % save
        cd(pData);
        titlename = sprintf('%s/%s composit %s',pData,figureprefix,gname{1});
        print(figure1,'-depsc', titlename); % save as (r600 will give better resolution)
        print(figure1,'-dtiff', titlename); % save as (r600 will give better resolution)
        close;
    end
end


%% script done

fprintf('\nDONE\n');






