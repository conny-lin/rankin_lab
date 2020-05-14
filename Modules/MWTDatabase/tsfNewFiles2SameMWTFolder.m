function [pdlist] = tsfNewFiles2SameMWTFolder(pDataSource,pDataBase,tsfonly)
% pDataBase = '/Volumes/COBOLT/MWT2';
% pDataSource = '/Volumes/COBOLT/Keep';




%% check 
addpath('/Users/connylin/Dropbox/MATLAB/Function_Library_Public');
% mwtnameSearch = '\<\d{8}[_]\d{6}\>)';
% excludeList = {'Analysis'; 'MatlabAnalysis'};

if nargin == 2
    tsfonly = 0;
end

%% get MWT data from checking folder
addpath('/Users/connylin/Dropbox/RL/Code/Modules/MWTDatabase');
[pMWTA,fMWTA] = getpMWT(pDataSource,1);

%% get MWT data from database folder
p = pDataBase;
addpath('/Users/connylin/Dropbox/RL/Code/Modules/MWTDatabase');
[pMWTD,fMWTD] = getpMWT(pDataBase,1);



%% transfer files
% pTemp = '/Volumes/COBOLT/MWT_ZipTemp';
startN = 1;
endN = numel(pMWTA);
reportJump = 100:100:numel(pMWTA);
fprintf('\n\nstarting from x=%d\n',startN);
fileExclusion = '(trinitySummary)|(rasterData)|(README[.]dat)|(Groups_1Set[.]mat)|(SwanLake)|(tapN)|(swanlake)';

for x = startN:endN
    %% A
    if sum(reportJump == x) > 0
        fprintf('Processing %d/%d MWT files\n',x,numel(pMWTA));
    end
    pmwta = pMWTA{x};
    fmwta = fMWTA{x};
    % see if anything in the folder
    [fa,pa] = dircontent(pmwta);
    if isempty(fa) == 1
        fprintf('[%s]: no analysis files\n',fmwta);
    else % if not empty, delete .*
        i = regexpcellout(fa,'\<[.]\w*');
        cellfun(@delete,pa(i))
        pa(i) = []; fa(i) = [];  
        
        %% D find coorsponding Data files
        i = ismember(fMWTD,fmwta);
        if sum(i) > 1
            error('more than one match');
        elseif sum(i) == 0 % if folder does not exist, transfer the whole folder
            fprintf('-no [%s] found in database, transfer folder\n',fmwta);
            ps = pmwta;
            pd = regexprep(pmwta,pDataSource,pDataBase);
            pd = fileparts(pd);
            if isdir(pd) == 0; mkdir(pd); end
            if tsfonly == 1
                movefile(pmwta,pd,'f')
            else
                [~,pslist] = dircontent(ps);
                pdlist = regexprep(pslist,pDataSource,pDataBase)
                pdlist{1}
                pslist{1}
                for f = 1:numel(pslist)
                    copyfile(pslist{f},pdlist{f})
                end
            end
        else % if MWT folder already exist, transfer files only
            pmwtd = pMWTD{i};
            % get data files
            [fd,pd] = dircontent(pmwtd);
            i = regexpcellout(fd,'\<[.]\w*');
            cellfun(@delete,pd(i))
            fd(i) = []; pd(i) = [];
            % find missing files
            i = ismember(fa,fd);
            fa(i) = []; pa(i) = [];
            if isempty(fa) == 0
                fprintf('-found new files from [%s] in source, transfer database\n',fmwta);
                ps = pa;
                pd = regexprep(ps,pmwta,pmwtd);
                if tsfonly == 1
                    fprintf('%d/%d. moving [%d] files: %s\n',x,endN, numel(ps),regexprep(pmwta,pDataSource,''));
                    cellfun(@movefile,ps,pd,cellfunexpr(ps,'f'));
                else
                    fprintf('%d/%d. copying [%d] files: %s\n',x,endN, numel(ps),regexprep(pmwta,pDataSource,''));
                    cellfun(@copyfile,ps,pd,cellfunexpr(ps,'f'));
                end
            end
            rmdir(pmwta,'s');
        end
        
    end
end

fprintf('DONE\n');

return
    

    %%
    
        % check consistency of files
        
        
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




%%             % see if different prefix
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


   




