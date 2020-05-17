function [Legend,pMWTpass,pMWTfailed] = chormaster5(option,pMWTS,varargin)
%% [Legend,pMWTcS,fval,chorscript] = chormaster4(option,pMWTS)

%% UPDATES
% r20151016
%     validateMWTcontent function
%     the same folder must contain .jar files required for javapath
% r20160220
%     added varargin to allow output of chorscripts
%
%

%% settings 
convert2mat = 0;
deletedat = 0;
displayopt = true;
legendonly = false;
pMWTpass = {};
pMWTfailed = {};

%% process varargin
vararginProcessor

%% CHOR OUTPUT LEGENDS ----------------------------------------------------
s = {'t time -- always the first column unless included again ';
'f frame -- the frame number ';
'p persistence -- length of time object is tracked ';
'D id -- the object ID ';
'n number -- the number of objects tracked ';
'N goodnumber -- the number of objects passing the criteria given ';
'e area -- body area';
'm midline -- length measured along the curve of object ';
'M morphwidth -- mean width of body about midline ';
'w width -- width of the rectangle framing the body';
'W relwidth -- instantaneous width/average width ';
'l length -- measured along major axis, not curve of object ';
'L rellength -- instantaneous length/average length ';
'a aspect -- length/width ';
'A relaspect -- instantaneous aspect/average aspect ';
'k kink -- head/tail angle difference from body (in degrees) ';
'c curve -- average angle (in degrees) between body split into 5 segments ';
's speed -- speed of movement';
'S angular -- angular speed ';
'b bias -- fractional excess of time spent moving one way ';
'P pathlen -- distance traveled forwards (backwards=negative) ';
'd dir -- consistency of direction of motion ';
'x loc_x -- x coordinate of object (mm) ';
'y loc_y -- y coordinate of object (mm) ';
'u vel_x -- x velocity (mm/sec) ';
'v vel_y -- y velocity (mm/sec) ';
'o orient -- orientation of body (degrees, only guaranteed modulo pi) ';
'r crab -- speed perpendicular to body orientation ';
'1 tap -- whether a tap (stimulus 1) has occurred ';
'2 puff -- whether a puff (stimulus 2) has occurred';
'3 stim3 -- whether the first custom stimulus has occurred ';
'4 stim4 -- whether the second custom stimulus has occurred. ';
'^ :max -- maximum value ';
'_ :min -- minimum value ';
'# :number -- number of items considered in this statistic ';
'- :median -- median value ';
'* :std -- standard deviation ';
':sem :sem -- standard error ';
':var :var -- variance ';
':? :exists -- 1 if the value exists, 0 otherwise ';
':p25 :p25 -- 25th percentile ';
':p75 :p75 -- 75th percentile ';
':jitter :jitter -- estimate of measurement precision '};


%% JAVA ARGUMENTS----------------------------------------------------------
% path to java programs
javapath = fileparts(mfilename('fullpath')); % get path of this function
b = blanks(1); % blank
% call java 
javacall = 'java -jar'; 
% RAM limit
javaRAM = '-Xmx8G'; 
javaRAM7G = '-Xmx7G';
% .jar paths
beethoven = ['''',javapath,'/Beethoven_v2.jar','''']; % call beethoven 
chor = ['''',javapath,'/Chore_1.3.0.r1035.jar','''']; % call chor 
% chor calls 
map = '--map';
% settings 
pixelsize = '-p 0.027'; 
speed = '-s 0.1'; 
mintime = '-t 20'; 
nall = '-N all';
minmove = '-M 2'; 
shape = '--shadowless -S';
% plugins 
preoutline = '--plugin Reoutline::exp';  
prespine = '--plugin Respine';

% plugins (reversals) 
revbeethoven_trv = '--plugin MeasureReversal::tap::dt=1::collect=0.5::postfix=trv';
revignortap_sprevs = '--plugin MeasureReversal::postfix=sprevs';
rev_ssr = '--plugin MeasureReversal::tap::collect=0.5::postfix=ssr';

% dat output collection
odrunkposture = '-O drunkposture -o nNslwakb';
odrunkposture2 = '-O drunkposture2 -o nNslwakbcemM';
oconny = '-O conny -o 1nee#e*ss#s*SS#S*ll#l*LL#L*ww#w*aa#a*mm#m*MM#M*kk#k*bb#b*pp#p*dd#d'; % Conny's 
obeethoven = '-o nNss*b12M'; % standard for Beethoven
oshanespark = '-O shanespark -o nNss*b12M'; % standard for Beethoven
oevan = '-O evan -o nNss*b12'; % Evan's dat output
oevanall = '-O evanall -N all -o nNss*b12';
oswanlakeall = '-O swanlakeall -N all -o tnNemMawlkcspbd1';
oswanlake = '-O swanlake -o tnNemMawlkcspbd1e#m#M#a#w#l#k#c#s#p#b#d#e-m-M-a-w-l-k-c-s-p-b-d-e*m*M*a*w*l*kvc*s*p*b*d*';
% Trinity app analysis for rastor plot and spontaneous locomotion
otrinity = '-O trinity -N all -o nNss*b12xyMmeSakcr'; 
ostarfish = '-O starfish -N all -o nNss*b12xyMmeSakcr';
ogangnam = '-O gangnam -N all -o DpmcobdPsSruvxy1';


%% MAKE JAVA SYNTAX (chorescript) -----------------------------------------
chorscript = {};
matval = {};
switch option
    case 'Gangnam'
       chorscript{1} = [javacall,b,javaRAM,b,chor,b,pixelsize,b,speed,b,...
            mintime,b,minmove,b,shape,b,ogangnam,b,...
            preoutline,b,prespine,b]; 
        fval = {'*gangnam*'}; 
        matval = {'Gangnam.mat'};
    case 'StarFish'
        chorscript{1} = [javacall,b,javaRAM,b,chor,b,pixelsize,b,speed,b,...
            mintime,b,minmove,b,shape,b,ostarfish,b,...
            preoutline,b,prespine,b]; 
        fval = {'*starfish*'};
    case 'Trinity'
        outputoption = otrinity;
        pluginoption = revbeethoven_trv;
        chorscript{1} = [javacall,b,javaRAM,b,chor,b,pixelsize,b,speed,b,...
            mintime,b,minmove,b,nall,b,shape,b,...
            outputoption,b,preoutline,b,prespine,b,pluginoption,b];        
        fval = {'*trinity.*.dat';'*.trv'};
        matval = {'trinitySummary.mat'};
    case 'TrinityOnly'
        outputoption = otrinity;
        chorscript{1} = [javacall,b,javaRAM,b,chor,b,pixelsize,b,speed,b,...
            mintime,b,minmove,b,nall,b,shape,b,...
            outputoption,b,preoutline,b,prespine,b];        
        fval = {'*trinity.*.dat'};
        matval = {'trinitySummary.mat'};
    case 'LadyGaGa'
        chorscript{1} = [javacall,b,javaRAM,b,chor,b,pixelsize,b,speed,b,...
            mintime,b,minmove,b,shape,b,oevan,b,preoutline,b,prespine,b,...
            revignortap_sprevs,b]; 
        fval = {'*evan.dat';'*.sprevs'};
    case 'DrunkPosture'
%         chorscript{1} = [javacall,b,javaRAM,b,chor,b,pixelsize,b,speed,b,...
%             mintime,b,minmove,b,shape,b,oevan,b,odrunkposture,b,....
%             preoutline,b,prespine,b]; 
%         fval = {'*drunkposture.dat'};
        chorscript{1} = [javacall,b,javaRAM,b,chor,b,pixelsize,b,speed,b,...
            mintime,b,minmove,b,shape,b,oevan,b,odrunkposture2,b,....
            preoutline,b,prespine,b]; 
        fval = {'*drunkposture2.dat'};        
    case 'ShaneSpark'
        chorscript{1} = [javacall,b,javaRAM,b,chor,b,pixelsize,b,speed,b,...
            mintime,b,minmove,b,shape,b,oshanespark,b,preoutline,b,...
            prespine,b,revbeethoven_trv,b]; 
        fval = {'*.trv';'*shanespark.dat'};
    case 'DrunkPostureOnly'
        chorscript{1} = [javacall,b,javaRAM,b,chor,b,pixelsize,b,speed,b,...
            mintime,b,minmove,b,shape,b,odrunkposture2,b,....
            preoutline,b,prespine,b]; 
        fval = {'*drunkposture2.dat'};  
    case 'Beethoven'
        chorscript{1} = [javacall,b,javaRAM,b,chor,b,pixelsize,b,speed,b,...
            mintime,b,minmove,b,shape,b,obeethoven,b,preoutline,b,...
            prespine,b,revbeethoven_trv,b]; 
        fval = {'*.trv'};
    case 'BeethovenOnly'
        chorscript{1} = [javacall,b,javaRAM,b,chor,b,pixelsize,b,speed,b,...
            mintime,b,minmove,b,shape,b,obeethoven,b,preoutline,b,...
            prespine,b,revbeethoven_trv,b]; 
        fval = {'*.trv'};
    case 'AnnaPavlova'
        chorscript{1} = [javacall,b,javaRAM,b,chor,b,pixelsize,b,speed,b,...
            mintime,b,minmove,b,shape,b,oshanespark,b,preoutline,b,...
            prespine,b,revbeethoven_trv,b]; 
        fval = {'*.trv';'*shanespark.dat'};
    case 'Rastor'
        chorscript{1} = [javacall,b,javaRAM,b,chor,b,pixelsize,b,speed,b,...
            mintime,b,minmove,b,shape,b,oevanall,b,preoutline,b,prespine,b];
        fval = {'*evanall*'};
    case 'SwanLake'
        chorscript{1} = [javacall,b,javaRAM,b,chor,b,pixelsize,b,speed,b,...
            mintime,b,minmove,b,shape,b,oswanlakeall,b,preoutline,b,...
            prespine,b,revbeethoven_trv,b];   
        chorscript{2} = [javacall,b,javaRAM,b,chor,b,pixelsize,b,speed,b,...
            mintime,b,minmove,b,shape,b,oswanlake,b,preoutline,b,...
            prespine,b]; 
        fval = {'*.trv';'*swanlake.dat'; '*swanlakeall*'};
    case 'Glee'
        chorscript{1} = [javacall,b,javaRAM,b,chor,b,pixelsize,b,speed,b,...
            mintime,b,minmove,b,shape,b,oswanlake,b,preoutline,b,...
            prespine,b,revbeethoven_trv,b];   
        fval = {'*.trv';'*swanlake.*.dat'};
    otherwise
        error('option does not exist in chore master');
end
% validate chorescript is cell
if iscell(chorscript) ==0; error 'chorescript must be in cell array'; end


% % get fval
% fval = {};
% for x = 1:numel(chorscript)
%    str = chorscript{x};
%    a = regexp(str,'(?<=[-]O\s{1})\w*(?=[-]N\s{1}all)','match');
%    for y = 1:numel(a)
%        fval = [favl;['*.',a{y},'.*.dat']];
%    end
%    a = regexp(str,'(?<=[-]O\s{1})\w*(?![-]N\s{1}all)','match');
%    for y = 1:numel(a)
%        fval = [fval;['*.',a{y},'.dat']];
%    end
% %    a = regexp(str,'(?<![-]O\s{1})[-]o');
% %    for y = 1:numel(a)
% %        fval = [fval;{'*.dat'}];
% %    end
%    a = regexp(str,'(?<=postfix[=])\w*','match');
%    if ~isempty(a); fval = [fval;['*.',char(a)]]; end
% end
% fval


%% MAKE LEGEND %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
T = table;
a = regexpcellout(s,' ','split');
T.short = a(:,1);
T.long = a(:,2);
a = regexpcellout(s,'(?<=(--)\s).*','match');
a = regexprep(a,'\s\>','');
a = regexprep(a,'[.]\>','');
T.description = a;
LegendTable = T;

% translate
Legend = struct;
for x = 1:numel(chorscript)
    script = chorscript{x};
    a = regexp(script,'(?<=(-)o\s{1})\S+','match');
    aa = regexp(script,'(?<=(-)O\s{1})\S+','match');
    for y = 1:numel(a)
       codes = a{y};
       c = cellstr(codes');
       [i,j] = ismember(c,LegendTable.short);
       b = cell(size(c));
       c(i) = LegendTable.long(j(i));
       k = find(regexpcellout(c,'[:]'));
       c(k-1) = cellfun(@strcat,c(k-1),c(k),'UniformOutput',0);
       c(k) = [];
       Legend.(aa{y}) = [{'time'};c];
    end
end

if legendonly
    return
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% CHECK CHOR HISTORY -----------------------------------------------------
if displayopt; fprintf('Checking existing chor outputs\n'); end
% find pMWT has mat file output already
if ~isempty(matval)
    [pMWT_Mat0,v1] = checkchoroutput(pMWTS,matval,'displayopt',displayopt);
    pMWT_Mat1 = pMWTS(~v1);
else
    pMWT_Mat1 = {};
    pMWT_Mat0 = pMWTS;
end
% find pMWT has text file output already
[pMWTcS,~] = checkchoroutput(pMWT_Mat0,fval,'displayopt',displayopt);


%% clean out . files ------------------------------------------------------
pF = pMWTcS;
recycle('on');
for x = 1:numel(pF)
    [fn,p] = dircontent(pF{x});
    i = regexpcellout(fn,'\<[.]');
    if sum(i)>0
        fprintf('deleting temp files starting with "."...');
        cellfun(@delete,p(i))
        fprintf('done\n')
    end
end


%% RUN CHOR ---------------------------------------------------------------
if isempty(pMWTcS) == 0
    fprintf('Chor %d MWT files: ',numel(pMWTcS));
    nFile = numel(pMWTcS);
    for x = 1:nFile
        % report status
        pD = pMWTcS{x};
        [~,mwtfn] = fileparts(pD); 
        if displayopt; fprintf('\nexamining %d/%d MWTfile [%s]...',x,nFile,mwtfn);end
        % run chor
        if displayopt; fprintf('chor-ing [%s]...\n',mwtfn); end
        
        file = strcat('''',pD,''''); 
        for y = 1:numel(chorscript) 
            if displayopt
                system([chorscript{y} file], '-echo'); 
            else
                system([chorscript{y} file]); 
            end
        end  
    end
    if displayopt; fprintf('\nChor Completed\n'); end
end


%% CONFIRM CHOR FILES -----------------------------------------------------
[pMWTfailed,~] = checkchoroutput(pMWTcS,fval,'displayopt',displayopt);
if displayopt
if isempty(pMWTfailed)== false
    fprintf('\nThese failed:\n');
    disp(char(unique(pMWTfailed)));
end
end
% get pass files
pMWTpass = pMWTS(~ismember(pMWTS,pMWTfailed));


%% CONVERT TO MAT FILES ---------------------------------------------------
if convert2mat
    switch option
       case 'Gangnam'
            convertchorNall2mat(pMWTpass,'gangnam','Gangnam');
       case 'Trinity'
            convertTrinityDat2Mat(pMWTpass,deletedat); 
       case 'TrinityOnly'
            convertTrinityDat2Mat(pMWTpass,deletedat); 
       otherwise
    end
end



%% FINISH
if displayopt; fprintf('\nChormaster5 DONE\n\n'); end

end


%% HELP FILES ---------------------------------------
% Options: 
%   -?  --help               This message (use -? output for help on output type) 
%       --body-length-units  Speeds are in units of body lengths (default is mm) 
%       --from               Time from which to read data (in seconds, default 0) 
%       --graph              Bring up GUI to graph population data 
%       --header             Write tab-delimited description of each data column 
%   -I (--interactive)       Bring up GUI (same as --graph) 
%       --in                 Only use data points inside specified shape 
%       --ignore-outside-triggers   Ignore all data except near explicit triggers 
%   -m (--minimum-move-mm)   How far an object must move (in mm) to count 
%   -M (--minimum-move-body)   (same thing, except unit is object-lengths) 
%       --minimum-biased     If object travels this far, it's mostly forwards 
%       --map                Use GUI to display the data as a browsable map 
%   -n (--id)                Only use listed object IDs (use commas: -n 1,5,22) 
%   -N (--each-id)           Write one output file for each ID listed 
%       --no-output          Don't write any output 
%       --no-repeat          Remove any frames that appear to be repeated 
%       --out                Data must be outside specified shape. 
%   -o (--output)            Write specified output data (-? output for syntax) 
%   -O (--output-name)       Add an identifier to output 
%   -p (--pixelsize)         Size of one pixel, in mm 
%       --plugin             Use plugin; --plugin help gives generic help 
%       --prefix             Specify data file prefix explicitly 
%   -q (--quiet)             Don't print progress information to console 
%   -s (--speed-window)      Time window (in seconds) to average velocity 
%   -S (--segment)           Shape analysis of path: lines, arcs, etc. 
%       --shadowless         Only count objects after they move a body length 
%       --skip-zeros         Omit timepoints with zero objects found 
%       --spine-from-outline (Re)compute spine more robustly given outline 
%   -t (--minimum-time)      How long an object must last (in seconds) to count 
%   -T (--output-rate)       Time between output data points (in seconds) 
%       --to                 Time after which to ignore data (in seconds) 
%       --target             Place all output in specified directory (must exist) 
%       --trigger            Report a stimulus-triggered average to .trig file 
%       --trig-only          Only write triggered averages, not regular output 
%       --who                Print out object ID numbers that pass criteria 
% Format: 
%   directory must contain a MWT .summary file 
%   A .zip file containing the data can be specified instead of the directory. 
%     The corresponding directory will be created for output purposes. 
%   -m,M,p,s,t,--from,--to expect a floating-point value as an argument 
%   -O name turns output from prefix.dat to prefix.name.dat 
%     If only one -O is given, it will change the .pos file name also. 
%     If multiple -O's are given, only .dat files are changed, and there must be 
%       the same number of -o's and -O's (and will correspond in order) 
%   --trigger is followed by the duration of the averaging window (in seconds), 
%     a comma, and then comma-separated list containing either the time at which 
%     to trigger or the tap, puff, stim3, or stim4 keywords followed by a colon, 
%     the time before to take a measurement, a colon, and the time after to 
%     take a measurement.  (Numbers may be left blank; colons are required.) 
%     Multiple trigger statements are okay (each adds more columns to the file). 
%   -n or --id can be entered multiple times, and/or can contain multiple id 
%     numbers; all IDs are accumulated.  Numbers must be separated by commas 
%     with no spaces.  IDs that do not exist or fail criteria are excluded. 
%     The -N or --each-id variant appends a five-digit object ID number to 
%     the prefix, and creates one set of files for each object. 
%     -N all means output separately every object meeting the criteria. 
%     -n and -N are not compatible.  Use only one. 
%   --in and --out should be followed by either a center and radius (circle) 
%     as x,y,r, or two corners of a rectangle as x1,y1,x2,y2. 
% Examples: 
%   --trigger 1.0,5,tap:0.25:0.5,750 will average from 5-6 s after 
%     the start of recording, from 1.25 to 0.25 s before each tap, from 0.5 s 
%     to 1.5 s after each tap, and once more from 750-751 s. 
%   --trigger 0.2,tap::0.2 will average from 0.2 to 0.4 seconds after each tap 
%   --trigger 0.5,tap:0: will average from 0.5s before to 0 s before each tap 
%   --in 1,1,100,50 --out 25,25,5 would only take data from an elongated 
%     rectangle with a hole missing from its left side. 
% crankin@Leviathan:~$ java -Xmx1500m -jar '/home/crankin/Desktop/Chore_1.3.0.r1035.jar' --help -? output 
% Choreography 1.3.0 build 1035 
% Usage:  java Choreography [options] directory 
%      or java -jar Chore.jar [options] directory 
% Format: 
%   -o requires an argument specifying columns (separate long form with commas) 
%   all -- same as ftnNpsSlLwWaAmkbcd1234 
    % time & frame
        % t time -- always the first column unless included again 
        % f frame -- the frame number 
        % p persistence -- length of time object is tracked 
    % object numer
        % D id -- the object ID 
        % n number -- the number of objects tracked 
        % N goodnumber -- the number of objects passing the criteria given 
    % object body size
        % e area 
        % m midline -- length measured along the curve of object 
        % M morphwidth -- mean width of body about midline 
    % posture
        % w width 
        % W relwidth -- instantaneous width/average width 
        % l length -- measured along major axis, not curve of object 
        % L rellength -- instantaneous length/average length 
        % a aspect -- length/width 
        % A relaspect -- instantaneous aspect/average aspect 
        % k kink -- head/tail angle difference from body (in degrees) 
        % c curve -- average angle (in degrees) between body split into 5 segments 
    % movement
        % s speed 
        % S angular -- angular speed 
        % b bias -- fractional excess of time spent moving one way 
        % P pathlen -- distance traveled forwards (backwards=negative) 
        % d dir -- consistency of direction of motion 
        % x loc_x -- x coordinate of object (mm) 
        % y loc_y -- y coordinate of object (mm) 
        % u vel_x -- x velocity (mm/sec) 
        % v vel_y -- y velocity (mm/sec) 
        % o orient -- orientation of body (degrees, only guaranteed modulo pi) 
        % r crab -- speed perpendicular to body orientation 
    % stimulus
        % 1 tap -- whether a tap (stimulus 1) has occurred 
        % 2 puff -- whether a puff (stimulus 2) has occurred
        % 3 stim3 -- whether the first custom stimulus has occurred 
        % 4 stim4 -- whether the second custom stimulus has occurred. 

%   The output items can be followed by the statistic to report 
%     (default is to output the mean) 
%     ^ :max -- maximum value 
%     _ :min -- minimum value 
%     # :number -- number of items considered in this statistic 
%     - :median -- median value 
%     * :std -- standard deviation 
%       :sem -- standard error 
%       :var -- variance 
%     ? :exists -- 1 if the value exists, 0 otherwise 
%       :p25 -- 25th percentile 
%       :p75 -- 75th percentile 
%       :jitter -- estimate of measurement precision 
%   Long format items need at least one comma (add trailing comma if needed) 
% Examples: 
%   -o a_ and -o aspect:min, are the same thing 
%   -o fnww* and -o frame,number,width,width:std also are the same 
%   -o xy will give positions (useful in conjunction with -N option) 
%   -o uv will give velocity vectors (also useful with -N) 
%   -o CCCCC will run through five different plugin-computed quantities, in order 
%     (advancing to the next plugin when the previous has computed all it can) 
%   -o CC*CC* will run through two but compute mean and SD of each. 
% 

%% MOVED CODE

