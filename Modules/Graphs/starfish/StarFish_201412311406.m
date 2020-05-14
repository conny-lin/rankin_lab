

%% define paths
% program path
pProgram = '/Users/connylin/OneDrive/MATLAB/Functions_Developer';
pJava = '/Users/connylin/OneDrive/MATLAB/Java_Programs';
% add program path
addpath(pProgram);
% data path
pData = '/Volumes/ParahippocampalGyrus/MWT_Data';
% output
pO = '/Users/connylin/OneDrive/Dance_Output';


%% choose a set of data to play
ExpName = '20141116B_JS_300s0x0s0s_24hrsPE_nahMT14480';
% ExpName = '20141106C_SM_100s30x10s10s';
pExp = [pData,'/',ExpName];


%% USER INPUTS
display ' ';
display 'Name your analysis output folder';
AnalysisName = input(': ','s');


%% GENERATE CHOR JAVA SYNTAX (chorescript)

% GENERATE JAVA ARGUMENTS
% path to java programs
%javapath = [strrep(userpath,pathsep,''),'/MATLAB MWT/SubFun_Java'];
% javapath = [pProgram,'/Java'];
javapath = pJava;
choreversion = 'Chore_1.3.0.r1035.jar';

b = blanks(1); % blank
% call java 
javacall = 'java -jar'; javaRAM = '-Xmx7G'; javaRAM7G = '-Xmx7G';
beethoven = ['''',javapath,'/Beethoven_v2.jar','''']; % call beethoven 
chor = ['''',javapath,'/',choreversion,'''']; % call chor 
% chor calls 
map = '--map';
% settings 
pixelsize = '-p 0.027'; speed = '-s 0.1'; 
mintime = '-t 20'; minmove = '-M 2'; shape = '--shadowless -S';
% plugins 
preoutline = '--plugin Reoutline::exp';  
prespine = '--plugin Respine';
% plugins (reversals) 
revbeethoven_trv = '--plugin MeasureReversal::tap::dt=1::collect=0.5::postfix=trv';
revignortap_sprevs = '--plugin MeasureReversal::postfix=sprevs';
rev_ssr = '--plugin MeasureReversal::tap::collect=0.5::postfix=ssr';

% dat output 
odrunkposture = '-O drunkposture -o nNslwakb';
odrunkposture2 = '-O drunkposture2 -o nNslwakbcemM';

oconny = '-o 1nee#e*ss#s*SS#S*ll#l*LL#L*ww#w*aa#a*mm#m*MM#M*kk#k*bb#b*pp#p*dd#d'; % Conny's 
obeethoven = '-o nNss*b12M'; % standard for Beethoven
oshanespark = '-O shanespark -o nNss*b12M'; % standard for Beethoven
oevan = '-O evan -o nNss*b12'; % Evan's dat output
oevanall = '-O evanall -N all -o nNss*b12';
oswanlakeall = '-O swanlakeall -N all -o tnNemMawlkcspbd1';
oswanlake = '-O swanlake -o tnNemMawlkcspbd1e#m#M#a#w#l#k#c#s#p#b#d#e-m-M-a-w-l-k-c-s-p-b-d-e*m*M*a*w*l*kvc*s*p*b*d*';
ostarfish = '-O starfish -o nNss*b12xyMmeSakcr';

chorscript{1} = [javacall,b,javaRAM,b,chor,b,pixelsize,b,speed,b,...
            mintime,b,minmove,b,'-N all',b,shape,b,ostarfish,b,...
            preoutline,b,prespine,b]; 
fval = {'*starfish*'};


%% Choreography
% get all MWT folder paths
[~,pMWT] = dircontentmwt(pExp);

% validate if need to chor
needChorVal = true(size(pMWT));
for x = 1:numel(pMWT);
    % this only checks the first fval
    if isempty(dircontent(pMWT{x},fval{1})) == 0 
        needChorVal(x) = false;
    end
end
pMWTc = pMWT(needChorVal);

% chor
if isempty(pMWTc) == 0
display(sprintf('Need to chor %d MWT files',numel(pMWTc)));
str = 'Chor-ing MWTfile [%s]...';
for x = 1:numel(pMWTc); 
    [~,fn] = fileparts(pMWTc{x}); file = strcat('''',pMWTc{x},''''); 
    display(sprintf(str,fn));
    for x = 1:numel(chorscript) 
        system([chorscript{x} file], '-echo'); 
    end  
end
end
display 'Chor Completed';  


%% create output folder
display ' ';
display 'Name your analysis output folder';
name = ['StarFish_',generatetimestamp,'_',AnalysisName];
% cd(pO);
pOFolder = [pO,'/',name];
if isdir(pOFolder) == 0; mkdir(pOFolder); end


%% STARFISH
% For Matlab set current folder to the time-date stamped experiment folder,
% then use this script (start and finish variables define the time window 
% of interest):

% DEFINE START-FINISH TIME SERIES
startTime = [40 [100:10:290]];
finishTime = [60 [110:10:300]];

% LOAD DATA TO MATLAB (EFFICIENCY REASON)
display ' ';
display 'converting raw data...';
for MWTfnX = 1:numel(pMWT)
    cd(pMWT{MWTfnX}); % cd to MWT folder
    [~,a] = fileparts(pMWT{MWTfnX});
    display(sprintf('converting [%s]',a));
    dirData = dir('*starfish*');
    for i = 1:numel(dirData);
        DataImport{MWTfnX,i} = dlmread(dirData(i).name);
    end
end

% save Data
% cd(pOFolder);
% save('matlab.mat','DataImport');


%% Making graphs
display ' ';
display 'making starfish graphs';
for t = 1:numel(startTime)
    
    for MWTfnX = 1:numel(pMWT)
        [p,MWTfn] = fileparts(pMWT{MWTfnX}); % get MWT file name
        
        % create output sub folder
        [~,GroupName] = fileparts(p);
        % make group folder if don't have one
        pOGroup = [pOFolder,'/',GroupName];
        if isdir(pOGroup) == 0; mkdir(pOGroup); end

        
        % evan's starfish code ---------------------
%         cd(pMWT{MWTfnX}); % cd to MWT folder

        displaceTog = [0];

        start = startTime(t);
        finish = finishTime(t);
%         dirData = dir('*starfish*');
        figure1 = figure('Visible','off');
        hold on

        k=1;
        
        for i = 1:sum((cellfun(@numel,DataImport(MWTfnX,:)) >0));

%         for i = 1:numel(dirData);
%             storedData = dlmread(dirData(i).name);
            storedData = DataImport{MWTfnX,i};
            
            if storedData(1,1) < start && storedData(end,1) > finish;
                k=1+k;
                irTimes = storedData(:,1) <start |storedData(:,1)>finish;
                storedData(irTimes,:)=[];

                conx = isnan(storedData(:,9));
                storedData(conx,:)=[];

                cony = isnan(storedData(:,10));
                storedData(cony,:)=[];

                storedData(:,11)= (storedData(1,9)-storedData(:,9));
                storedData(:,12)= (storedData(1,10)-storedData(:,10));

                randCol = [rand rand rand];

                plot(storedData(:,11),storedData(:,12),'color',randCol,'LineWidth',1.5);

                displace = sqrt(storedData(end,11)*storedData(end,11) + storedData(end,12)*storedData(end,12));

                displaceTog = [displaceTog;displace];

            end
            clear conx cony irTimes storedData displace
        end
    % evan's starfish code ---------------------

    % save figure
    figureName = [MWTfn,'[',num2str(startTime(t)),'_', num2str(finishTime(t)),']'];
    savefigepsOnly(figureName,pOGroup);
    end

end

display 'Starfish done';














