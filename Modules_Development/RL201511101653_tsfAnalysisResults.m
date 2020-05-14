
%% check 
addpath('/Users/connylin/Dropbox/MATLAB/Function_Library_Public');
mwtnameSearch = '\<\d{8}[_]\d{6}[.]zip\>';
pCheck = '/Volumes/COBOLT/MWT_dup_to_delete';
excludeList = {'Analysis'; 'MatlabAnalysis'};



%% get zipped data
pHome = '/Volumes/COBOLT/MWT';
[~,~,~,pEf] = dircontent(pHome);
% get group folder
[~,~,fG,pG] = cellfun(@dircontent,pEf,'UniformOutput',0);
fG = celltakeout(fG);
pG = celltakeout(pG);
[~,~,fMWT,pMWT] = cellfun(@dircontent,pG,'UniformOutput',0);
pMWT = celltakeout(pMWT,'multirow');
fMWT = celltakeout(fMWT,'multirow');
i = regexpcellout(fMWT,'(\<\d{8}[_]\d{6}\>)|(\<\d{8}[_]\d{6}[.]zip\>)');
pMWT(~i) = [];
fMWT(~i) = [];
fMWTD = fMWT;
pMWTD = pMWT;



%% get MWT data from checking folder
pA = '/Users/connylin/Dropbox/Lab/MWT_Analysis';
[fE,pE,fEf,pEf] = dircontent(pA);
% get group folder
[~,~,fG,pG] = cellfun(@dircontent,pEf,'UniformOutput',0);
fG = celltakeout(fG);
pG = celltakeout(pG);
[~,~,fMWT,pMWT] = cellfun(@dircontent,pG,'UniformOutput',0);
pMWT = celltakeout(pMWT,'multirow');
fMWT = celltakeout(fMWT,'multirow');
i = regexpcellout(fMWT,'\<\d{8}[_]\d{6}\>');
pMWT(~i) = [];
fMWT(~i) = [];
fMWTA = fMWT;
pMWTA = pMWT;


%% transfer files
% pTemp = '/Volumes/COBOLT/MWT_ZipTemp';
startN = 1600;
endN = numel(pMWTA);
reportJump = 100:100:numel(pMWTA);
fprintf('\n\nstarting from x=%d\n',startN);
fileExclusion = '(trinitySummary)|(rasterData)|(README[.]dat)|(Groups_1Set[.]mat)|(SwanLake)|(tapN)|(swanlake)';
for x = startN:endN
    if sum(reportJump == x) > 0
        fprintf('Processing %d/%d MWT files\n',x,numel(pMWTA));
    end
    pmwta = pMWTA{x};
    fmwta = fMWTA{x};
    
    %% find coorsponding Data files
    i = ismember(fMWTD,fmwta);
    if sum(i) > 1
        error('more than one match or no match');
    end
    if sum(i) == 1;
        pmwtd = pMWTD{i};
    
        %% check consistency of files
        [fa,pa] = dircontent(pmwta);
        if isempty(fa) == 1
            fprintf('[%s]: no analysis files\n',fmwta);
        else
            i = regexpcellout(fa,'\<[.]\w*');
            cellfun(@delete,pa(i))
            pa(i) = [];
            fa(i) = [];    
%             % see if different prefix
%             f = fa;
%             ftest = f; % eliminate some particular files
%             ftest(regexpcellout(ftest,fileExclusion)) = [];
% 
%             prefix = regexprep(ftest(regexpcellout(ftest,'([.]trv)')),'([.]trv)','');
%             if isempty(prefix) == 1
%                 prefix = regexprep(ftest(regexpcellout(ftest,'([.]png)')),'([.]png)','');
%             end
%             if isempty(prefix) == 1
%                 prefix = regexprep(ftest(regexpcellout(ftest,'([.]set)')),'([.]set)','');
%             end
%             if isempty(prefix) == 1
%                 a = regexp(ftest,'[.]');
%                 b = nan(size(a,1),1);
%                 for ai = 1:size(a,1)
%                     b(ai,1) = a{ai}(1);
%                 end
%                 if sum(diff(b)) == 0;
%                     prefix = f{1}(1:b(1)-1);
%                     if sum(regexpcellout(ftest,prefix)) ~= numel(ftest)
%                         error('stop, too many prefix');
%                     end
%                 else
%                     f
%                     error('dot at the wrong spot');
%                 end
%             elseif numel(prefix) > 1
%                 error('stop, too many prefix for a');
%             end
%             prefixA = prefix;
%             ftestA = ftest;
% 
%             i = regexpcellout(ftest,prefix);
%             if sum(i) ~= numel(ftest); 
%                 warning('more than one prefix')
%                 fprintf('std prefix of a: %s\n',char(prefix));
%                 disp(char(f(~i)));
%                 opt = input('override? y-1 no-0: ');
%                 if opt == 0;         
%                     startN = x+1;
%                     error('stop x = %d',x); 
%                 end
%             end

            %% get data files
            [fd,pd] = dircontent(pmwtd);
            i = regexpcellout(fd,'\<[.]\w*');
            cellfun(@delete,pd(i))
            fd(i) = [];
            pd(i) = [];


%             % see if different prefix
%             f = fd;
%             ftest = f; % eliminate some particular files
%             ftest(regexpcellout(ftest,fileExclusion)) = [];
% 
%             prefix = regexprep(dircontent(pmwtd,'*.set'),'[.]set','');
%             if isempty(prefix) == 1
%                 prefix = regexprep(dircontent(pmwtd,'*.summary'),'[.]summary','');
%             elseif numel(prefix) > 1
%                 error('stop, too many prefix')
%             end
%             prefixD = prefix;
%             ftestD = ftest;
%              i = regexpcellout(ftest,prefix);
%             if sum(i) ~= numel(ftest); 
%                 warning('more than one prefix')
%                 fprintf('std prefix of a: %s\n',char(prefix));
%                 disp(char(f(~i)));
%                 opt = input('override? y-1 no-0: ');
%                 if opt == 0;         
%                     startN = x+1;
%                     error('stop x = %d',x); 
%                 end
%             end
% 
%             %% check d and a prefix match
%             if strcmp(prefixA,prefixD) == 0
%                 error('analysis and data prefix do not match\n-a=%s\n-d=%s\n',char(prefixA),char(prefixD));
%             end

            %% find file not in data folder
            i = ~ismember(fa,fd);
            if sum(i) > 0;
                pD = regexprep(pa(i),pA,pHome);
                pS = pa(i);
                fprintf('%d. copying [%d] files: %s\n',x, numel(pD),regexprep(pmwta,pA,''));
                cellfun(@copyfile,pS,pD,cellfunexpr(pD,'f'));
                startN = x+1;
            end
        end

    end
        
end

return

%% checks
p1 = '/Users/connylin/Dropbox/Lab/MWT_Analysis/20120120C_CL_100s30x10s10s_LiquidTsf/N2_TsfLiquid/20120120_180214';

a = dir([p1,'/N2_5x3_t96h20C_100s30x10s10s_C0116ae.set'])

% p2 = '/Volumes/COBOLT/MWT/20120120C_CL_100s30x10s10s_LiquidTsf/N2_TsfLiquidM9/20120120_163322';
% b = dir([p2,'/N2_5x3_t96h20C_100s30x10s10s_C0116ba.set'])

%%
pMWTD(ismember(fMWTD,'20120120_170216'))
pMWTA(ismember(fMWTA,'20120120_170216'))




