
%% check 
addpath('/Users/connylin/Dropbox/MATLAB/Function_Library_Public');
pDeleteArchive = '/Volumes/COBOLT/MWT_dup_to_delete';

excludeList = {'Analysis'; 'MatlabAnalysis'};
mwtnameSearch = '\<\d{8}[_]\d{6}\>';


%% find MWT files found in analysis but not in data
p = '/Users/connylin/Dropbox/Lab/MWT_Analysis';
[~,~,fE,pE] = dircontent(p);
pEA = pE;
fEA = fE;
[~,~,~,pG] = cellfun(@dircontent,pE,'UniformOutput',0);
pG = celltakeout(pG,'multirow');
[~,~,fMWT,pMWT] = cellfun(@dircontent,pG,'UniformOutput',0);
pMWT = celltakeout(pMWT,'multirow');
fMWT = celltakeout(fMWT,'multirow');
i = regexpcellout(fMWT,mwtnameSearch);
fMWT(~i) = [];
pMWT(~i) = [];
% load to A
pMWTA = pMWT;
fMWTA = fMWT;
fMWTAU = unique(fMWTA);


p = '/Volumes/COBOLT/MWT';
[~,~,fE,pE] = dircontent(p);
pED = pE;
fED = fE;
[~,~,~,pG] = cellfun(@dircontent,pE,'UniformOutput',0);
pG = celltakeout(pG,'multirow');
[~,~,fMWT,pMWT] = cellfun(@dircontent,pG,'UniformOutput',0);
pMWT = celltakeout(pMWT,'multirow');
fMWT = celltakeout(fMWT,'multirow');
i = regexpcellout(fMWT,mwtnameSearch);
fMWT(~i) = [];
pMWT(~i) = [];
% load to A
pMWTD = pMWT;
fMWTD = fMWT;
fMWTDU = unique(fMWTD);

%% check unique
if numel(fMWTDU) ~= numel(fMWTD)
    error('duplicate data MWT');
end
if numel(fMWTAU) ~= numel(fMWTA)
    error('duplicate analysis MWT');
end

%% remove prefix
p = '/Users/connylin/Dropbox/Lab/MWT_Analysis';
pMWTATest = regexprep(pMWTA,p,'');
p = '/Volumes/COBOLT/MWT';
pMWTDTest = regexprep(pMWTD,p,'');

%% get problem plates
i = ismember(pMWTDTest,pMWTATest);
pMWTDProblem = pMWTDTest(~i);
size(pMWTDProblem)

%% create table for B
p = pMWTDProblem;
[pG,fMWT] = cellfun(@fileparts,p,'UniformOutput',0);
[pE,fG] = cellfun(@fileparts,pG,'UniformOutput',0);
[~,fE] = cellfun(@fileparts,pE,'UniformOutput',0);
T = table;
T.pMWT = p;
T.fMWT = fMWT;
T.fG = fG;
T.fE = fE;
TD =T;

% create table for A
p = pMWTATest;
[pG,fMWT] = cellfun(@fileparts,p,'UniformOutput',0);
[pE,fG] = cellfun(@fileparts,pG,'UniformOutput',0);
[~,fE] = cellfun(@fileparts,pE,'UniformOutput',0);
T = table;
T.pMWT = p;
T.fMWT = fMWT;
T.fG = fG;
T.fE = fE;
TA =T;


%% put aside D MWT without reference in A

i = ~ismember(TD.fMWT,TA.fMWT);
TD_noRef = TD(i,:);
TD(i,:) = [];

%% keep only Analysis reference with TD mwt
i = ~ismember(TA.fMWT,TD.fMWT);
TA(i,:) = [];

% sort
TD = sortrows(TD,{'fMWT'});
TA = sortrows(TA,{'fMWT'});




%% check if difference fG
TD.fGA = TA.fG;
TD.fGmatch = true(size(TD,1),1);
for x = 1:size(TD,1)
    TD.fGmatch(x) = strcmp(TD.fG{x},TD.fGA{x});
end
% 6 folders have right group folder match
sum(TD.fGmatch) 

%% check Exp folder
TD.fEA = TA.fE;
%%
TD.fEmatch = true(size(TD,1),1);
for x = 1:size(TD,1)
    TD.fEmatch(x) = strcmp(TD.fE{x},TD.fEA{x});
end
sum(~TD.fEmatch) 
% 15 has wrong exp folder match

%% manully fix exp 
TD(~TD.fEmatch,:) = [];

%%
TD.fGmatch = true(size(TD,1),1);
for x = 1:size(TD,1)
    TD.fGmatch(x) = strcmp(TD.fG{x},TD.fGA{x});
end
% no folders have right group folder match
sum(TD.fGmatch) 

%% change folders
% group into folders
fGDU = unique(TD.fG);

%%
TD.pGD = cellfun(@fileparts,TD.pMWT,'UniformOutput',0);
pGDU = unique(TD.pGD);
pD = '/Volumes/COBOLT/MWT';
%%


for g = 270:numel(pGDU)
    psG = pGDU{g};
    i = ismember(TD.pGD,psG);
    fsG = unique(TD.fG(i));
    fdG = unique(TD.fGA(i));
    fsE = unique(TD.fE(i));
    
    if numel(fdG) > 1 || numel(fsG) > 1
        T2 = TD(i,:);
        for g2 = 1:numel(fdG)
            fdG2 = fdG{g2};
            i = ismember(T2.fGA,fdG2);
            T3 = T2(i,:);
            psF = cellfun(@strcat,cellfunexpr(T3.pMWT,pD),T3.pMWT,'UniformOutput',0);
            psG = [pD,char(unique(T3.pGD))];
            pd = [fileparts(psG),'/',char(fdG2)];
            if isdir(pd) == 0; mkdir(pd); end
            cellfun(@movefile,psF,cellfunexpr(psF,pd),cellfunexpr(psF,'f'));
        end
        % remove old folder
        [a] = dircontent(psG);   
        if isempty(a) == 1;
            rmdir(psG,'s');
        else
            disp(psG)
            error('folder not empty g = %d',g)
        end
        
        
    else
        fprintf('from [%s]\n--moving [%s] to [%s]\n',char(fsE),char(fsG),char(fdG))
%         a = input('ok?(y-1,n-0): ');
%         if a ~=1, warning('manual stop'); return; end
        
        % move the whole folder
        psG = [pD,TD.pGD{i}];
        pd = [fileparts(psG),'/',char(fdG)];
        if isdir(pd) == 0; mkdir(pd); end
        [~,psF] = dircontent(psG);
        cellfun(@movefile,psF,cellfunexpr(psF,pd),cellfunexpr(psF,'f'));
       % remove old folder
       [a] = dircontent(psG);   
       if isempty(a) == 1;
           rmdir(psG,'s');
       else
           disp(psG)
           error('folder not empty g = %d',g)
       end
        
        
    end
end
 

%%

for x = 5:size(TD,1)
    
    ps = [pD,TD.pMWT{x}];
    pd = fileparts([pD,TA.pMWT{x}]);
    if isdir(pd) == 0; mkdir(pd); end
    movefile(ps,pd,'f');
    % check if empty
    psF = fileparts(ps);
    c = dircontent(psF);
    if isempty(c) == 1; 
        rmdir(psF); 
    else 
    end
end




%%
i = ismember(TD.fE,'20141026C_JS_3600s0x0s0s_acute200_RB1025');
TD(i,:)
%%
TD(i,:) = [];




%% match
pD = '/Volumes/COBOLT/MWT';

for x = 5:size(TD,1)
    
    ps = [pD,TD.pMWT{x}];
    pd = fileparts([pD,TA.pMWT{x}]);
    if isdir(pd) == 0; mkdir(pd); end
    movefile(ps,pd,'f');
    % check if empty
    psF = fileparts(ps);
    c = dircontent(psF);
    if isempty(c) == 1; 
        rmdir(psF); 
    else 
    end
end


%%
unique(TD_noRef.fE)



%% take out problem E that's not in Analysis folder
p = pMWTA;
[pG,fMWT] = cellfun(@fileparts,p,'UniformOutput',0);
[pE,fG] = cellfun(@fileparts,pG,'UniformOutput',0);
[~,fE] = cellfun(@fileparts,pE,'UniformOutput',0);

fEPU(~ismember(fEPU,fE)) = [];

%% take out MWT files not have analysis folder reference
[pG,fMWT] = cellfun(@fileparts,pMWTDProblem,'UniformOutput',0);
[pE,fG] = cellfun(@fileparts,pG,'UniformOutput',0);
[~,fE] = cellfun(@fileparts,pE,'UniformOutput',0);
i = ismember(fE,fEPU);
pMWTDPCheck = pMWTDProblem(i);

p = pMWTDPCheck;
[pG,fMWT] = cellfun(@fileparts,p,'UniformOutput',0);
[pE,fG] = cellfun(@fileparts,pG,'UniformOutput',0);
[~,fE] = cellfun(@fileparts,pE,'UniformOutput',0);






%% unify group name

startN = 1;
%%
for x = startN:numel(pED)
    [~,~,fGD,pGD] = dircontent(pED{x});
    fGD(ismember(fGD,excludeList)) = [];

    i = ismember(fEA,fED{x});
    if sum(i) ~= 1
        disp(char(fGD))
        startN = x+1;
        error('x=%d, no exp match in analysis: %s\n',x,fED{x});
    else
        [~,~,fGA,pGA] = dircontent(pEA{i});
        fGA(ismember(fGA,excludeList)) = [];
        if sum(ismember(fGD,fGA)) ~= numel(fGD)
            disp(char(fGA))
            disp(char(fGD))
            
            
            %%
            if analysis has less, ignore
                
            ismember(
            
            %%
            startN = x+1;
            error('x=%d group name not match\n-Exp:%s\n',x,fED{x});
        else
            fprintf('x=%d group name match, check MWT folders\n-Exp:%s\n',x,fED{x});

        end
    end
    
end






