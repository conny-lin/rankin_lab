function syncDatabases(pMatch,pHome)
%% check 
addpath('/Users/connylin/Dropbox/Code/Matlab/Library/General');
excludeList = {'Analysis'; 'MatlabAnalysis'};
% pMatch = '/Volumes/COBOLT/MWT_Data_Archive_Zip'; 
pMatchType = 'zip';
% pHome = '/Volumes/COBOLT/MWT';
pHomeType = 'folder';

%% find MWT files found in analysis but not in data
[~,~,fE,pE] = dircontent(pMatch);
pEA = pE;
fEA = fE;
[~,~,~,pG] = cellfun(@dircontent,pE,'UniformOutput',0);
pG = celltakeout(pG,'multirow');
switch pMatchType
    case 'folder'
        [~,~,fMWT,pMWT] = cellfun(@dircontent,pG,'UniformOutput',0);
        mwtnameSearch = '\<\d{8}[_]\d{6}\>';
    case 'zip'
        [fMWT,pMWT] = cellfun(@dircontent,pG,'UniformOutput',0);
        mwtnameSearch = '\<\d{8}[_]\d{6}[.]zip\>';

end
pMWT = celltakeout(pMWT,'multirow');
fMWT = celltakeout(fMWT,'multirow');
i = regexpcellout(fMWT,mwtnameSearch);
fMWT(~i) = [];
pMWT(~i) = [];
% load to A
pMWTA = pMWT;
fMWTA = fMWT;
fMWTAU = unique(fMWTA);


[~,~,fE,pE] = dircontent(pHome);
pED = pE;
fED = fE;
[~,~,~,pG] = cellfun(@dircontent,pE,'UniformOutput',0);
pG = celltakeout(pG,'multirow');
switch pHomeType
    case 'folder'
        [~,~,fMWT,pMWT] = cellfun(@dircontent,pG,'UniformOutput',0);
        mwtnameSearch = '\<\d{8}[_]\d{6}\>';
    case 'zip'
        [fMWT,pMWT] = cellfun(@dircontent,pG,'UniformOutput',0);
        mwtnameSearch = '\<\d{8}[_]\d{6}[.]zip\>';
end
pMWT = celltakeout(pMWT,'multirow');
fMWT = celltakeout(fMWT,'multirow');
i = regexpcellout(fMWT,mwtnameSearch);
fMWT(~i) = [];
pMWT(~i) = [];
% load to A
pMWTD = pMWT;
fMWTD = fMWT;
fMWTDU = unique(fMWTD);



%% check if all mwt foldres are unique
if numel(fMWTDU) ~= numel(fMWTD)
    error('duplicate data MWT');
end
if numel(fMWTAU) ~= numel(fMWTA)
    error('duplicate analysis MWT');
end



%% remove prefix
pMWTATest = regexprep(pMWTA,pMatch,'');
pMWTDTest = regexprep(pMWTD,pHome,'');
if strcmp(pMatchType,'zip') == 1
    pMWTATest = regexprep(pMWTATest,'[.]zip','');
end
if strcmp(pHomeType,'zip') == 1
    pMWTDTest = regexprep(pMWTDTest,'[.]zip','');
end


%% get problem plates
i = ismember(pMWTATest,pMWTDTest);
pMWTAProblem = pMWTATest(~i);
fprintf('\n%d/%d MWT files need to change names\n',size(pMWTAProblem,1),numel(pMWTATest))
if strcmp(pMatchType,'zip') == 1
    pMWTAProblem = cellfun(@strcat,pMWTAProblem,cellfunexpr(pMWTAProblem,'.zip'),'UniformOutput',0);
end


%% move files
for x= 1:numel(pMWTAProblem)
    pS = [pMatch,pMWTAProblem{x}];
    disp(pMWTAProblem{x});
    [pSG,fmwts] = fileparts(pS);
    i = ismember(fMWTD,fmwts);
    if sum(i) > 1; error('more than one mwt'); end
    
    if sum(i) == 1
        pD = regexprep(pMWTD{i},pHome,pMatch);
        disp(['to ',char(regexprep(pD,pMatch,''))]);
        pDG = fileparts(pD);
        if isdir(pDG) == 0; mkdir(pDG); end
        ps = pS;
        pd = pDG;
        movefile(ps,pd,'f');
        if strcmp(pMatchType,'folder')
            pGold = fileparts(pS);
            [~,pf,~,pF] = dircontent(pGold);
            i = regexpcellout(pF,'(\<\d{8}[_]\d{6}\>)|(\<\d{8}[_]\d{6}[.]zip\>)'); % find files with proper MWT names
            j = regexpcellout(pf,'(\<\d{8}[_]\d{6}\>)|(\<\d{8}[_]\d{6}[.]zip\>)'); % find files with proper MWT names
            if isempty(i) == 1 && isempty(pf) > 0; % if no more MWT files, move all files
                disp(ps); 
                cellfun(@movefile,pf,cellfunexpr(pf,pd),cellfunexpr(pf,'f'));
            end
        end
        if isempty(dircontent(pGold)) == 1; rmdir(pGold,'s'); end

    end
end





end
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% %% create table for B
% p = pMWTDProblem;
% [pG,fMWT] = cellfun(@fileparts,p,'UniformOutput',0);
% [pE,fG] = cellfun(@fileparts,pG,'UniformOutput',0);
% [~,fE] = cellfun(@fileparts,pE,'UniformOutput',0);
% T = table;
% T.pMWT = p;
% T.fMWT = fMWT;
% T.fG = fG;
% T.fE = fE;
% TD =T;
% 
% % create table for A
% p = pMWTATest;
% [pG,fMWT] = cellfun(@fileparts,p,'UniformOutput',0);
% [pE,fG] = cellfun(@fileparts,pG,'UniformOutput',0);
% [~,fE] = cellfun(@fileparts,pE,'UniformOutput',0);
% T = table;
% T.pMWT = p;
% T.fMWT = fMWT;
% T.fG = fG;
% T.fE = fE;
% TA =T;
% 
% 
% %% put aside D MWT without reference in A
% 
% i = ~ismember(TD.fMWT,TA.fMWT);
% TD_noRef = TD(i,:);
% TD(i,:) = [];
% 
% %% keep only Analysis reference with TD mwt
% i = ~ismember(TA.fMWT,TD.fMWT);
% TA(i,:) = [];
% 
% % sort
% TD = sortrows(TD,{'fMWT'});
% TA = sortrows(TA,{'fMWT'});
% 
% 
% 
% 
% %% check if difference fG
% TD.fGA = TA.fG;
% TD.fGmatch = true(size(TD,1),1);
% for x = 1:size(TD,1)
%     TD.fGmatch(x) = strcmp(TD.fG{x},TD.fGA{x});
% end
% % 6 folders have right group folder match
% sum(TD.fGmatch) 
% 
% %% check Exp folder
% TD.fEA = TA.fE;
% %%
% TD.fEmatch = true(size(TD,1),1);
% for x = 1:size(TD,1)
%     TD.fEmatch(x) = strcmp(TD.fE{x},TD.fEA{x});
% end
% sum(~TD.fEmatch) 
% % 15 has wrong exp folder match
% 
% %% manully fix exp 
% TD(~TD.fEmatch,:) = [];
% 
% %%
% TD.fGmatch = true(size(TD,1),1);
% for x = 1:size(TD,1)
%     TD.fGmatch(x) = strcmp(TD.fG{x},TD.fGA{x});
% end
% % no folders have right group folder match
% sum(TD.fGmatch) 
% 
% %% change folders
% % group into folders
% fGDU = unique(TD.fG);
% 
% %%
% TD.pGD = cellfun(@fileparts,TD.pMWT,'UniformOutput',0);
% pGDU = unique(TD.pGD);
% pD = '/Volumes/COBOLT/MWT';
% %%
% 
% 
% for g = 1:numel(pGDU)
%     psG = pGDU{g};
%     i = ismember(TD.pGD,psG);
%     fsG = unique(TD.fG(i));
%     fdG = unique(TD.fGA(i));
%     fsE = unique(TD.fE(i));
%     
%     if numel(fdG) > 1 || numel(fsG) > 1
%         T2 = TD(i,:);
%         for g2 = 1:numel(fdG)
%             fdG2 = fdG{g2};
%             i = ismember(T2.fGA,fdG2);
%             T3 = T2(i,:);
%             psF = cellfun(@strcat,cellfunexpr(T3.pMWT,pD),T3.pMWT,'UniformOutput',0);
%             psG = [pD,char(unique(T3.pGD))];
%             pd = [fileparts(psG),'/',char(fdG2)];
%             if isdir(pd) == 0; mkdir(pd); end
%             cellfun(@movefile,psF,cellfunexpr(psF,pd),cellfunexpr(psF,'f'));
%         end
%         % remove old folder
%         [a] = dircontent(psG);   
%         if isempty(a) == 1;
%             rmdir(psG,'s');
%         else
%             disp(psG)
%             error('folder not empty g = %d',g)
%         end
%         
%         
%     else
%         fprintf('from [%s]\n--moving [%s] to [%s]\n',char(fsE),char(fsG),char(fdG))
% %         a = input('ok?(y-1,n-0): ');
% %         if a ~=1, warning('manual stop'); return; end
%         
%         % move the whole folder
%         psG = [pD,TD.pGD{i}];
%         pd = [fileparts(psG),'/',char(fdG)];
%         if isdir(pd) == 0; mkdir(pd); end
%         [~,psF] = dircontent(psG);
%         cellfun(@movefile,psF,cellfunexpr(psF,pd),cellfunexpr(psF,'f'));
%        % remove old folder
%        [a] = dircontent(psG);   
%        if isempty(a) == 1;
%            rmdir(psG,'s');
%        else
%            disp(psG)
%            error('folder not empty g = %d',g)
%        end
%         
%         
%     end
% end
%  
% 
% %%
% 
% for x = 5:size(TD,1)
%     
%     ps = [pD,TD.pMWT{x}];
%     pd = fileparts([pD,TA.pMWT{x}]);
%     if isdir(pd) == 0; mkdir(pd); end
%     movefile(ps,pd,'f');
%     % check if empty
%     psF = fileparts(ps);
%     c = dircontent(psF);
%     if isempty(c) == 1; 
%         rmdir(psF); 
%     else 
%     end
% end
% 
% 
% 
% 
% %%
% i = ismember(TD.fE,'20141026C_JS_3600s0x0s0s_acute200_RB1025');
% TD(i,:)
% %%
% TD(i,:) = [];
% 
% 
% 
% 
% %% match
% pD = '/Volumes/COBOLT/MWT';
% 
% for x = 5:size(TD,1)
%     
%     ps = [pD,TD.pMWT{x}];
%     pd = fileparts([pD,TA.pMWT{x}]);
%     if isdir(pd) == 0; mkdir(pd); end
%     movefile(ps,pd,'f');
%     % check if empty
%     psF = fileparts(ps);
%     c = dircontent(psF);
%     if isempty(c) == 1; 
%         rmdir(psF); 
%     else 
%     end
% end
% 
% 
% %%
% unique(TD_noRef.fE)
% 
% 
% 
% %% take out problem E that's not in Analysis folder
% p = pMWTA;
% [pG,fMWT] = cellfun(@fileparts,p,'UniformOutput',0);
% [pE,fG] = cellfun(@fileparts,pG,'UniformOutput',0);
% [~,fE] = cellfun(@fileparts,pE,'UniformOutput',0);
% 
% fEPU(~ismember(fEPU,fE)) = [];
% 
% %% take out MWT files not have analysis folder reference
% [pG,fMWT] = cellfun(@fileparts,pMWTDProblem,'UniformOutput',0);
% [pE,fG] = cellfun(@fileparts,pG,'UniformOutput',0);
% [~,fE] = cellfun(@fileparts,pE,'UniformOutput',0);
% i = ismember(fE,fEPU);
% pMWTDPCheck = pMWTDProblem(i);
% 
% p = pMWTDPCheck;
% [pG,fMWT] = cellfun(@fileparts,p,'UniformOutput',0);
% [pE,fG] = cellfun(@fileparts,pG,'UniformOutput',0);
% [~,fE] = cellfun(@fileparts,pE,'UniformOutput',0);
% 
% 
% 
% 
% 
% 
% %% unify group name
% 
% startN = 1;
% %%
% for x = startN:numel(pED)
%     [~,~,fGD,pGD] = dircontent(pED{x});
%     fGD(ismember(fGD,excludeList)) = [];
% 
%     i = ismember(fEA,fED{x});
%     if sum(i) ~= 1
%         disp(char(fGD))
%         startN = x+1;
%         error('x=%d, no exp match in analysis: %s\n',x,fED{x});
%     else
%         [~,~,fGA,pGA] = dircontent(pEA{i});
%         fGA(ismember(fGA,excludeList)) = [];
%         if sum(ismember(fGD,fGA)) ~= numel(fGD)
%             disp(char(fGA))
%             disp(char(fGD))
%             
%             
%             %%
%             if analysis has less, ignore
%                 
%             ismember(
%             
%             %%
%             startN = x+1;
%             error('x=%d group name not match\n-Exp:%s\n',x,fED{x});
%         else
%             fprintf('x=%d group name match, check MWT folders\n-Exp:%s\n',x,fED{x});
% 
%         end
%     end
%     
% end






