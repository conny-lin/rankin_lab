function [Stats_byPlates,Data] = importTrv(pMWT)
%% .TRV OUTPUT LEGENDS
% output legends
% time
% # worms already moving backwards (can't score) 
% # worms that didn't reverse in response 
% # worms that did reverse in response 
% mean reversal distance (among those that reversed) 
% standard deviation of reversal distances 
% standard error of the mean 
% minimum 
% 25th percentile 
% median 
% 75th percentile 
% maximum 
% mean duration of reversal (also among those that reversed 
% standard deviation of duration 
% standard error of the mean 
% minimum 
% 25th percentile 
% median 
% 75th percentile 
% maximum


%% IMPORT .TRV % revised 20151009
A = pMWT;
for m = 1:size(A,1);
    [~,p] = dircontent(A{m,1},'*.trv'); 
    % if there is no .trv
    if isempty(p) == 1
        A{m,2} = {};     
    else       
        % validate trv output format
        pt = p{1};
        fileID = fopen(pt,'r');
        d = textscan(fileID, '%s', 2-1, 'Delimiter', '', 'WhiteSpace', '');
        fclose(fileID);
        % read trv
        if strcmp(d{1}{1}(1),'#') ==1 % if trv file is made by Beethoven
            a = dlmread(pt,' ',5,0); 
        else % if trv file is made by Dance
            a = dlmread(pt,' ',0,0);
        end
        A{m,2} = a(:,[1,3:5,8:10,12:16,19:21,23:27]); % index to none zeros
    end
end
Data = A;

% MWTfnImport = A;
% MWTSet.Data.Import = MWTfnImport;

% legend
% L = {'time','N?','N_NoResponse','N_Reversed','?','RevDist'    };

%% CHECK TAP CONSISTENCY
[r,c] = cellfun(@size,Data(:,2),'UniformOutput',0);
rn = celltakeout(r,'singlenumber');
rfreq = tabulate(rn);
rcommon = rfreq(rfreq(:,2) == max(rfreq(:,2)),1);
str = 'Common tap number = %d';
display(sprintf(str,rcommon));
rproblem = rn ~= rcommon;

if sum(rproblem)~=0;
    % sort data
    Data_tapProblem = Data(rproblem,:);
    Data = Data(~rproblem,:);
    
    % get paths and name of problem files
    pMWTfP = Data_tapProblem(:,1);
    [~,MWTfnP] = cellfun(@fileparts,pMWTfP,'UniformOutput',0);

    str = 'The following MWT did not have the same tap(=%d)';
    display(sprintf(str,rcommon)); 
    disp(MWTfnP);
    
%     % export report
%     T = cell2table(MWTfnP,'VariableNames',{'MWT'});
%     cd(pSave);
%     writetable(T,'MWT_plate_with_missing_taps.csv');
    
    display 'Removing from analysis...'; 
%     MWTSet.RawBad = Data(rproblem,:);
%     Data = Data(~rproblem,:);
%     MWTfnOK = MWTfn(~rproblem);
%     pMWTfOK = pMWTf(~rproblem);    
% 
%     % reconstruct
% %     [MWTSet.MWTfG] = reconstructMWTfG(pMWTf);
%     MWTSet.MWTInfo.pMWT = pMWTfOK;
%     MWTSet.MWTInfo.MWTfn = MWTfnOK;
%     MWTSet.MWTInfo.pMWTBadTap = pMWTf(rproblem); 
%     [~,g] = cellfun(@fileparts,...
%         cellfun(@fileparts,pMWTfOK,'UniformOutput',0),...
%         'UniformOutput',0);
%     MWTSet.MWTInfo.GroupName = g;

end

%% MAKING SENSE OF TRV 
% get data
% D = A; 
A = Data;
B = struct;
[~,n] = cellfun(@fileparts,A(:,1),'UniformOutput',0);
% B.MWTfn = MWTSet.MWTInfo.MWTfn;
B.MWTfn = n;
B.pMWT = A(:,1);
B.Raw = A(:,2);
% get group names
p = B.pMWT;
p = cellfun(@fileparts,p,'UniformOutput',0);
[p,n] = cellfun(@fileparts,p,'UniformOutput',0);
B.gname = n;
[~,n] = cellfun(@fileparts,p,'UniformOutput',0);
B.expname = n;

% pMWTf = MWTSet.MWTInfo.pMWT;
pMWTf = A(:,1);

% indexes of .trv
ind.RevDur = 13;
ind.RevDist = 5;
% calculation

for m = 1:size(pMWTf,1);
    % X = tap time
    % B.X.TapTime(:,m) = A{m,2}(:,1);
    B.X(:,m) = A{m,2}(:,1);   
    % basic caluations
    B.N.NoResponse(:,m) = A{m,2}(:,3);
    B.N.Reversed(:,m) = A{m,2}(:,4);  
    B.N.TotalN(:,m) = B.N.Reversed(:,m)+B.N.NoResponse(:,m);

    %% N
    n = B.N.TotalN(:,m);
    N = B.N.TotalN(:,m);
    N(n < 1) = NaN;

    % Frequency
    B.Y.RevFreq(:,m) = B.N.Reversed(:,m)./N;
    % variance can not be calculated at this point
    B.E.RevFreq(:,m) = NaN(size(B.N.Reversed(:,m))); %  can only be zero
    B.SD.RevFreq(:,m) = NaN(size(B.N.Reversed(:,m)));
    % B.Y.RevFreq(:,m) = B.N.Reversed(:,m)./B.N.TotalN(:,m);
    % B.Y.RevFreqStd(:,m) = B.Y.RevFreq(:,m)/B.Y.RevFreq(1,m);

    % Reversal Duration
    B.Y.RevDur(:,m) = A{m,2}(:,ind.RevDur);
    B.E.RevDur(:,m) = A{m,2}(:,ind.RevDur+1)./B.N.Reversed(:,m);
    B.SD.RevDur(:,m) = A{m,2}(:,ind.RevDur+1)./B.N.Reversed(:,m);

    % Reversal Speed = RevDist/RevDur
    % Distance [disabled]
    RevDist(:,m) = A{m,2}(:,ind.RevDist); 
%     B.SD.RevDist(:,m) = A{m,2}(:,ind.RevDist+1);
%     B.E.RevDist(:,m) = A{m,2}(:,ind.RevDist+1)./B.N.Reversed(:,m);
    % B.Y.RevDistStd(:,m) = B.Y.RevDist(:,m)/B.Y.RevDist(1,m);
    % B.Y.SumRevDist(:,m) = B.Y.RevDist(:,m).*B.N.Reversed(:,m);    
    B.Y.RevSpeed(:,m) = RevDist(:,m)./B.Y.RevDur(:,m); 
    B.E.RevSpeed(:,m) = NaN(size(B.Y.RevSpeed(:,m))); 
    B.SD.RevSpeed(:,m) = NaN(size(B.Y.RevSpeed(:,m)));
end



Stats_byPlates = B;
% MWTSet.Data.Raw = Raw;