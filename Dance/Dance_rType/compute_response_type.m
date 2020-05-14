function [T1, RTclean, RT, RTN,legend] = compute_response_type(Response, Baseline, varargin)

%% default
pSave = '';
vararginProcessor

%% get baseline direction
baseline_dir = Baseline(:,end);
baseline_dir(baseline_dir > 0) = 1;
baseline_dir(baseline_dir < 0) = -1;


%% get rid of baseine data in diff direction from the last baseline dir
B = Baseline;
B(B(baseline_dir == 1,:) < 0) = NaN;
B(B(baseline_dir == -1,:) > 0) = NaN;
B(isnan(baseline_dir),:) = NaN;


%% calculate baseline range
baseline_max = max(B')';
baseline_maxM = repmat(baseline_max, 1, size(Response,2));
baseline_min = min(B')'; 
baseline_minM = repmat(baseline_min, 1, size(Response,2));


%% calculate response difference
dR = nan(size(Response));
dR(Response > baseline_maxM) = 1;
dR(Response < baseline_maxM) = -1;
dR(Response <= baseline_maxM & Response >= baseline_minM) = 0;


%% make decision
baseline_dirM = repmat(baseline_dir, 1, size(Response, 2));

% text
RT = cell(size(Response)); % declare output array
RT(dR == 0) = {'baseline'};
RT(dR == -1 & baseline_dirM == 1 & Response > 0) = {'decelerate forward'};
RT(dR == -1 & baseline_dirM == 1 & Response < 0) = {'accelerate forward'};
RT(dR == -1 & baseline_dirM == 1 & Response == 0) = {'pause'};
RT(dR == -1 & baseline_dirM == -1 & Response < 0) = {'accelerate reverse'};
RT(dR == -1 & baseline_dirM >= 0 & Response < 0) = {'reversal'};
RT(dR == 1 & baseline_dirM >= 0 & Response > 0) = {'accelerate forward'};
RT(dR == 1 & baseline_dirM < 0 & Response < 0) = {'decelerate reverse'};
RT(dR == 1 & baseline_dirM >= 0 & Response == 0) = {'pause'};


%% numeric legend
legend = table;
legend.id = [1:7]';
legend.name = {'accelerate forward';
    'accelerate reverse';
    'pause';
    'baseline';
    'decelerate forward';
    'decelerate reverse';
    'reverse'};

%% numeric
RTN = nan(size(RT));
RTN(dR == 0) = 4;
RTN(dR == -1 & baseline_dirM == 1 & Response > 0) = 5;
RTN(dR == -1 & baseline_dirM == 1 & Response < 0) = 1;
RTN(dR == -1 & baseline_dirM == 1 & Response == 0) = 3;
RTN(dR == -1 & baseline_dirM == -1 & Response < 0) = 2;
RTN(dR == -1 & baseline_dirM >= 0 & Response < 0) = 7;
RTN(dR == 1 & baseline_dirM >= 0 & Response > 0) = 1;
RTN(dR == 1 & baseline_dirM < 0 & Response < 0) = 6;
RTN(dR == 1 & baseline_dirM >= 0 & Response == 0) = 3;



%% condense response switch text
% find change data
a = diff(RTN')';

b = false(size(RTN,1), size(RTN,2)-1);
b(a ~=0 & ~isnan(a)) = true;

c = true(size(RTN,1),1);
c(isnan(RTN(:,1))) = false;

b = [c b];
validata = b;

% delete repeat response type
A = RT;
A(~validata) = {''};
RTclean = A;

A = strjoinrows(A,'-');
% clean 
A = regexprep(A,'[-]{2,}','-');
A = regexprep(A,'([-]\>)|(\<[-])','');
A = regexprep(A,'[-]', ' > ');


%% save if prompted
T = tabulate(A);

T1 = table;
T1.response = T(:,1);
T1.n = cell2mat(T(:,2));
T1.p = cell2mat(T(:,3));



















































