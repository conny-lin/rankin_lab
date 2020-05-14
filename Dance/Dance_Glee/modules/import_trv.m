function Data = import_trv(pMWT)
%% Data = import_trv(pMWT) (revised 20151126)
% import .trv chor output into table
% mwtpath = pMWT list
% trv = trv import: if no trv present, empty entry
% .trv can be made by beethoven or dance. the format will be different


%% legend
L = {'time',...
'N_alreadyRev',...% # worms already moving backwards (can''t score)',
'N_ForwardOrPause',...%# worms that didn''t reverse in response', 
'N_Rev',...%# worms that did reverse in response', 
'RevDis',...%mean reversal distance (among those that reversed)',
'RevDis_SD',...%standard deviation of reversal distances', 
'RevDis_SE',...%standard error of the mean', 
'RevDis_min',...%minimum',
'RevDis_25percentile',...%25th percentile',
'RevDis_median',... 
'RevDis_75percentile',...
'RevDis_max', ...
'RevDur',...%mean duration of reversal (also among those that reversed',
'RevDur_SD',...%standard deviation of duration', 
'RevDur_SE',...%standard error of the mean', 
'RevDur_min',...
'RevDur_25percentile',...
'RevDur_median', ...
'RevDur_75percentile',...
'RevDur_max'
};

% A = cell(size(pMWTS,1),2); A(:,1) = pMWTS; % potentially junk code
B = table;
B.mwtpath = pMWT;
B.data = cell(size(pMWT,1),1);
for m = 1:size(pMWT,1);
    [~,p] = dircontent(pMWT{m},'*.trv'); 
    % delete temperarary .*.trv files
    if size(p,1) ~= 1
        cellfun(@delete,p(regexpcellout(p,'\<[.]\w{1,}[.]trv\>')))
        [~,p] = dircontent(pMWT{m},'*.trv');
    end
    if isempty(p) == 0    
        % validate trv output format
        pt = p{1};
        fileID = fopen(pt,'r');
        d = textscan(fileID,'%s', 2-1,'Delimiter', '', 'WhiteSpace', '');
        fclose(fileID);
        % read trv
        if strcmp(d{1}{1}(1),'#') ==1 % if trv file is made by Beethoven
            a = dlmread(pt,' ',5,0); 
        else % if trv file is made by Dance
            a = dlmread(pt,' ',0,0);
        end
%         A{m,2} = a(:,[1,3:5,8:10,12:16,19:21,23:27]); % index to none zeros
        % convert to table
        b = array2table(a(:,[1,3:5,8:10,12:16,19:21,23:27]),'VariableNames',L);
        B.data{m} = b;
    end
end
% MWTfnImport = A;
% MWTSet.Data.Import = B;
Data = B;
