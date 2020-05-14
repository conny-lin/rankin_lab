function D = trinityExtractData_v1707(TrinityData,taps,beforetap,aftertap,var,varargin)

%% DEFAULTS 
frameInt = 0.2; % probably the lowest one can go
saveopt = 0;
displayopt = false;
% var = {'wormid','time','tap','speed','bias'};
D = {};

%% VARARGIN PROCESSOR
vararginProcessor;

%% get legend
D = extract_trinity2table(TrinityData,var);

%% Deal with time 
% find recording time
recordTimeu = unique(D.time);
frameintmax = max(diff(recordTimeu));
tapTimeu = unique(D.time(D.tap==1));

if numel(tapTimeu) ~= taps
   error('tap number in the file does not equal to taps requested'); 
end

isi = unique(round(diff(tapTimeu)));
if numel(isi)>1
    disp(isi)
    warning('this code can not accomodate variable ISI');
    return
end
assayStartTimes = tapTimeu-beforetap;
assayEndTimes = tapTimeu+aftertap;


% create frame numbers
assayTimesCentreAdj = -beforetap:frameInt:aftertap;



% record assay time period number
recordTimePeriodTble = table;
recordTimePeriodTble.time = recordTimeu;
% put in tap times
a = zeros(size(recordTimePeriodTble,1),1);
for ti = 1:numel(tapTimeu)
    i = recordTimePeriodTble.time == tapTimeu(ti);
    a(i) = ti;
end
recordTimePeriodTble.centretime = a; % center time is usually a tap


% enter period time
tval = nan(size(recordTimeu));
recordTimePeriodTble.assayperiod = nan(size(recordTimePeriodTble,1),1);
recordTimePeriodTble.frame = recordTimePeriodTble.assayperiod;
for ti = 1:numel(assayStartTimes)
    
    t1 = assayStartTimes(ti);
    tf = assayEndTimes(ti);
    
    recordi = recordTimeu >= t1 & recordTimeu <= tf; % get times wtihin assay window
    
    tval(recordi) = ti;
    
    T = recordTimePeriodTble(recordi,:);
    centre_time = T.time(T.centretime>0);
    framet = assayTimesCentreAdj+ centre_time;
    frametind = nan(size(T,1),1);
    for fi = 1:numel(framet)
        ta = framet(fi);
        [m,i] = min(abs(T.time-ta));
        frametind(i) = fi;
        
    end

    recordTimePeriodTble.frame(recordi) = frametind;    
end
recordTimePeriodTble.assayperiod = tval;
recordTimePeriodTble(isnan(recordTimePeriodTble.frame),:) = [];


%% remove data outside of assay time
D(~ismember(D.time,recordTimePeriodTble.time),:) = [];
D = innerjoin(D,recordTimePeriodTble);



end

