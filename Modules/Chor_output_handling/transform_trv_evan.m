function B = transform_trv_evan(Data)
%% transform_trv_evan(A)
% A is a cell array containing trv imports


%% MAKING SENSE OF TRV (r20151126)
% get data
A = Data.data;
pMWT = Data.mwtpath;
B = struct;
[~,n] = cellfun(@fileparts,pMWT,'UniformOutput',0);
B.MWTfn = n;
B.pMWT = pMWT;

% calculation
for m = 1:size(B.pMWT,1);
    % tap time
    B.X(:,m) = A{m}.time;   
    % basic caluations
    B.N.NoResponse(:,m) = A{m}.N_ForwardOrPause;
    B.N.Reversed(:,m) = A{m}.N_Rev;  
    B.N.TotalN(:,m) = A{m}.N_ForwardOrPause + A{m}.N_Rev;

    % N
    n = B.N.TotalN(:,m);
    N = B.N.TotalN(:,m);
    N(n < 1) = NaN;

    % Frequency
    B.Y.RevFreq(:,m) = B.N.Reversed(:,m)./N;
    % variance can not be calculated at this point
    B.E.RevFreq(:,m) = NaN(size(B.N.Reversed(:,m))); %  can only be zero
    B.SD.RevFreq(:,m) = NaN(size(B.N.Reversed(:,m)));

    % Reversal Duration
    B.Y.RevDur(:,m) = A{m}.RevDur;
    B.E.RevDur(:,m) = A{m}.RevDur_SE; 
    B.SD.RevDur(:,m) = A{m}.RevDur_SD;

    % Reversal Speed = RevDist/RevDur  
    B.Y.RevSpeed(:,m) = A{m}.RevDis./A{m}.RevDur; 
    B.E.RevSpeed(:,m) = NaN(size(B.Y.RevSpeed(:,m))); 
    B.SD.RevSpeed(:,m) = NaN(size(B.Y.RevSpeed(:,m)));
end
