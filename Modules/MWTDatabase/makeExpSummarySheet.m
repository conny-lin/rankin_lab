function makeExpSummarySheet(databasefile,pSave,varargin)
%% make experiment summary sheet

%% varargin
nStart = 1;
expname = {};
vararginProcessor;


addpath('/Users/connylin/Dropbox/Code/Matlab/Library/Graphs');
addpath('/Users/connylin/Dropbox/Code/Matlab/Library/General');

% pSave = '/Users/connylin/Dropbox/RL/Publication InPrep/PhD Dissertation/Chapters/2-Materials and Methods General/MWT experiments';

%% EXAMINE DATABASE
load(databasefile);
D = MWTDB.text;

% max MWT per exp (63)
a = tabulate(D.expname);
nmax_mwtperexp = max(cell2mat(a(:,2)));
fprintf('\nmax mwt per exp: %d\n',nmax_mwtperexp);

% max number of plate per group (9)
enu = unique(D.expname);
nmax_mwtpergroup = 4;
nmax_groupperexp = 4;
en_max = {};
gn_max_mwt = {};
for ei = 1:numel(enu)
   en = enu{ei};
   d = D(ismember(D.expname, enu{ei}),:);
   a = tabulate(d.groupname);
   if size(a,1) > nmax_groupperexp
       nmax_groupperexp = size(a,1);
       en_max = en;
   end
   n = max(cell2mat(a(:,2)));
   if n > nmax_mwtpergroup
       nmax_mwtpergroup = n;
   end
end
fprintf('max mwt per group: %d\n',nmax_mwtpergroup);
fprintf('max group per exp: %d\n',nmax_groupperexp);


%% print exp sheets
% bad files: 26, 27
ht = 11;
w = 8.5;
if isempty(expname)
    enu = unique(D.expname);
else
    enu = expname;
end

%%

for ei = nStart:numel(enu)
    close;
    en = enu{ei};
    fprintf('experiment %d/%d: %s...',ei,numel(enu),en)
    
    % get exp database
    d = D(ismember(D.expname,en),:);
    
    %% set up page
    fig1 = figure('Visible','off','Color',[1 1 1]);
    axes1 = axes('Parent',fig1,'ZColor',[1 1 1],'YColor',[1 1 1],'XColor',[1 1 1]); 
    hold on;
    fillPage(fig1);
    xlim([0 w]);
    ylim([0 ht]);
    %% experiment title
    en = enu{ei};
    x = w/2; y = ht; 
    txt = regexprep(en,'_',' ');
    text(x,y,txt,'FontSize',16,'HorizontalAlignment','center')
    % draw divider
    x = 0:w; y=y-0.2; yline = repmat(y,numel(x),1);
    line(x,yline,'Color','k');
    
    %% list strains used
    y = y-0.1;
    genotype = d.genotype;
    genotype(cellfun(@isempty,genotype)) = [];
    genotypeu = unique(genotype);
    if isempty(genotypeu) == 1
        txt = 'no genotype info found in strain database';
        text(0,y,txt,'FontSize',10,'HorizontalAlignment','left');
        y = y-0.2;

    else
        for gi = 1:numel(genotypeu)
            s = unique(d.strain(ismember(d.genotype,genotypeu{gi})));
            txt = sprintf('%s = %s',char(s),genotypeu{gi});
            text(0,y,txt,'FontSize',10,'HorizontalAlignment','left')        
            y = y-0.2;
        end        
    end

    % draw divider
    x = 0:w; yline = repmat(y,numel(x),1);
    line(x,yline,'Color','k');

    % list groups and number of MWT
    a = tabulate(d.groupname);
    x = 0; 
    for gi = 1:size(a,1)
        y = y-0.2;
        gn = regexprep(a{gi,1},'_',' ');
        n = a{gi,2};
        txt = sprintf('(%d) %s',n,gn);
        text(x,y,txt,'FontSize',10,'HorizontalAlignment','left')        
    end
    
    % draw divider
    x = 0:w; y = y-0.2; yline = repmat(y,numel(x),1);
    line(x,yline,'Color','k');
    
    %% get MWT information
    % sort by group
    y=y-0.2;
    d = sortrows(d,{'groupname'});
    mwtpath = d.mwtpath;
    for mwti = 1:numel(mwtpath)
        %% set y
        %% get mwt
        pmwt = mwtpath{mwti};
        [p,mwtn] = fileparts(pmwt);
        mwtn = regexprep(mwtn,'_','-');
        [~,gn] = fileparts(p);
        gn = regexprep(gn,'_',' ');
        %% # of blobs
        nblobs = numel(dircontent(pmwt,'*.blob*'));
        
        %% get info from .set file
        [fn,p] = dircontent(pmwt,'*.set');
        p(regexpcellout(fn,'\<[.]')) = [];
        rc = '?';
        if isempty(p) == 0
            delimiter = '\t';
            startRow = 1;
            endRow = inf;
            formatSpec = '%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%[^\n\r]';
            fileID = fopen(p{1},'r');
            dset = textscan(fileID, formatSpec, endRow(1)-startRow(1)+1, 'Delimiter', delimiter, 'HeaderLines', startRow(1)-1, 'ReturnOnError', false);
            for block=2:length(startRow)
                frewind(fileID);
                dataArrayBlock = textscan(fileID, formatSpec, endRow(block)-startRow(block)+1, 'Delimiter', delimiter, 'HeaderLines', startRow(block)-1, 'ReturnOnError', false);
                for col=1:length(dset)
                    dset{col} = [dset{col};dataArrayBlock{col}];
                end
            end
            fclose(fileID);
            
            a = cell(numel(dset{1,1}),numel(dset));
            for x = 1:numel(dset)
                a(:,x) = dset{x};
            end
            %% set info
            original_filepath = a{1,1}; % .set file path
            mwtsetname = a{1,2}; % .set file name
            % .set file runcond
            recordtime = str2double(a{1,21}); % total recording time
            preplate = str2double(a{1,23}); % preplate
            isi = str2double(a{1,24}); % isi
            tapN = str2double(a{1,25}); % tapN
            postrecord = recordtime-((isi*(tapN-1))+preplate); % posttap
            firsttap = str2num(a{end,1})/1000; % first tap time
            rc = sprintf('%ds%dx%ds%ds',preplate,tapN,isi,postrecord);
        end
        
        
        %% load summary file
        tapN_actual = NaN;
        % worm N
        n_max = NaN;
        n_min = NaN;
        n_median = NaN;        
        % worm size
        size_max = NaN;
        size_min = NaN;
        size_median = NaN;
        [fn,p] = dircontent(pmwt,'*.summary');
        p(regexpcellout(fn,'\<[.]')) = [];
        if isempty(p) == 0
            delimiter = ' ';
            formatSpec = '%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%[^\n\r]';
            fileID = fopen(p{1},'r');
            dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'MultipleDelimsAsOne', true,  'ReturnOnError', false);
            fclose(fileID);
            a = cell(numel(dataArray{1,1}),numel(dataArray));
            for x = 1:numel(dataArray)
                b = dataArray{x};
                a(1:numel(b),x) = b;
            end
            clear dataArray;
            
            %% get info from .summary file
            % tap number actual
            b = a(:,17);
            i = regexpcellout(b,'0x1');
            tapN_actual = sum(i);
            %% worm N
            b = a(:,4);
            b(cellfun(@isempty,b)) = {'0'};
            n = cellfun(@str2num,b);
            n_max = max(n);
            n_min = min(n);
            n_median = median(n);        
            %% worm size
            b = a(:,15);
            b(cellfun(@isempty,b)) = {'0'};
            s = cellfun(@str2num,b);
            s = cellfun(@str2num,b);
            size_max = max(s);
            size_min = min(s);
            size_median = median(s); 
        end
        
        %% text
%         mwtsetname = regexprep(mwtsetname,'_','-');
        txt = sprintf('[%s]  %s: blobs(%d) rc:%s(tap%d) N:%d %d %d. size:%.0f %.0f %.0f',...
            gn,mwtn,nblobs,rc,tapN_actual,n_min,n_median,n_max,size_min,...
            size_median,size_max);
        text(0,y,txt,'FontSize',10,'HorizontalAlignment','left')  
%         original_filepath = regexprep(original_filepath,'\','-');
%         txt = sprintf('    fn: %s path: %s',mwtsetname,original_filepath);
%         y = y-0.2;
%         text(0,y,txt,'FontSize',8,'HorizontalAlignment','left')  

        %% move y
        y = y-0.2;
    end
    
    cd(pSave);
    print(fig1,'-dpdf',en);
    fprintf('-- printed\n');
    close all;
    
end
































