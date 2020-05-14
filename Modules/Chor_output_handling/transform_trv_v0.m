function B = transform_trv_v0(Data)
%% transform_trv_evan(A)
% A is a cell array containing trv imports


%% MAKING SENSE OF TRV (r20151126)
% get data
A = Data.data;
pMWT = Data.mwtpath;
B = struct;
[~,n] = cellfun(@fileparts,pMWT,'UniformOutput',0);
% B.MWTfn = n;
% B.pMWT = pMWT;
B.RevFreq.pMWT = pMWT';
B.RevDur.pMWT = pMWT';
B.RevSpeed.pMWT = pMWT';
B.RevFreq.mwtname = n';
B.RevDur.mwtname = n';
B.RevSpeed.mwtname = n';
% calculation
for m = 1:size(pMWT,1);
    nRow = size(A{m}.N_Rev,1);
    emptyArray = nan(nRow,1);
    % Frequency
    B.RevFreq.N(:,m) = A{m}.N_Rev + A{m}.N_ForwardOrPause;
    B.RevFreq.time(:,m) = A{m}.time;
    B.RevFreq.mean(:,m) = A{m}.N_Rev ./ (A{m}.N_ForwardOrPause + A{m}.N_Rev);
    % variance can not be calculated at this point
    B.RevFreq.SE(:,m) = emptyArray; %  can only be zero
    B.RevFreq.SD(:,m) = emptyArray;

    % Reversal Duration
    B.RevDur.N(:,m) = A{m}.N_Rev;
    B.RevDur.time(:,m) = A{m}.time;
    B.RevDur.mean(:,m) = A{m}.RevDur;
    B.RevDur.SE(:,m) = A{m}.RevDur_SE; 
    B.RevDur.SD(:,m) = A{m}.RevDur_SD;

    % Reversal Speed = RevDist/RevDur 
    B.RevSpeed.N(:,m) = A{m}.N_Rev;
    B.RevSpeed.time(:,m) = A{m}.time;
    B.RevSpeed.mean(:,m) = A{m}.RevDis./A{m}.RevDur; 
    B.RevSpeed.SE(:,m) = emptyArray; 
    B.RevSpeed.SD(:,m) = emptyArray;
end
