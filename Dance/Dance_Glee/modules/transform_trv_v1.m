function B = transform_trv_v1(Data,varargin)
%% transform_trv_evan(A)
% A is a cell array containing trv imports

%% DEFAULTS AND VARARGIN
clean = 1;
% process varargin
vararginProcessor; 

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
    d = A{m};

    nRow = size(d.N_Rev,1);
    emptyArray = nan(nRow,1);
    
    % Frequency
    B.RevFreq.N(:,m) = d.N_Rev + d.N_ForwardOrPause;
    B.RevFreq.time(:,m) = d.time;
    B.RevFreq.mean(:,m) = d.N_Rev ./ (d.N_ForwardOrPause + d.N_Rev);
    % variance can not be calculated at this point
    B.RevFreq.SE(:,m) = emptyArray; %  can only be zero
    B.RevFreq.SD(:,m) = emptyArray;

    
    % Reversal Duration
    B.RevDur.N(:,m) = d.N_Rev;
    B.RevDur.time(:,m) = d.time;
    B.RevDur.mean(:,m) = d.RevDur;
    B.RevDur.SE(:,m) = d.RevDur_SE; 
    B.RevDur.SD(:,m) = d.RevDur_SD;
    

    % Reversal Speed = RevDist/RevDur 
    B.RevSpeed.N(:,m) = d.N_Rev;
    B.RevSpeed.time(:,m) = d.time;
    B.RevSpeed.mean(:,m) = d.RevDis./d.RevDur; 
    B.RevSpeed.SE(:,m) = emptyArray; 
    B.RevSpeed.SD(:,m) = emptyArray;
end

%% if cleaning speed is yes, then remove nan speed
if clean == 1
D = B;
msrlist = fieldnames(D);
for msri = 1:numel(msrlist)
    msr = msrlist{msri};
    d = D.(msr).mean;
    % remove nan data from dataset
    i = any(isnan(d));
    if sum(i) > 1
        iremove = i;
        a = D.(msr);
        fn = fieldnames(a);
        A = struct;
        for fni = 1:numel(fn)
            A.(fn{fni}) = a.(fn{fni})(:,~iremove);
        end
        D.(msr) = A; % put it back
    end
end
end
