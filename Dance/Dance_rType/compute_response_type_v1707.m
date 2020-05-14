function [rDecision,rDecisionText,legend] = compute_response_type_v1707(Baseline_speed,Baseline_bias, Response_speed, Response_bias, varargin)

%% 
%% default
vararginProcessor

%% baseline
% flip baseline
BS = flip(Baseline_speed,2);
BSBias = flip(Baseline_bias,2);
[BSC,bDir] = remove_data_after_dir_change(BSBias,BS);
% calculate baseline range
BSV = BSC.*bDir;
BSV(bDir==0,:) = BSC(bDir==0,:);
bVmax = nanmax(BSV')';
bVmin = nanmin(BSV')';

%% response
[RSC,rDir] = remove_data_after_dir_change(Response_bias,Response_speed);
% calculate velocity with pause speed included a positive values
m = nanmean(RSC')';
rV = m.*rDir;
rV(rDir==0) = m(rDir==0); % for pauses, still record the speed as positive


%% numeric legend
legend = table;
legend.id = [1:8]';
legend.name = {...
'pause';       
'no_response'  
'reverse';  
'reverse_acc'; 
'reverse_dec'; 
'forward_acc'; 
'forward';
'forward_dec'};    

%% make decision
A = nan(size(Response_speed,1),1);
% baseline rev, response rev, response rev faster than baseline rev
A(bDir==-1 & rDir==-1 & rV < bVmin) = legend.id(ismember(legend.name,'reverse_acc')); 
% baseline rev, response rev, response rev slower than baseline rev
A(bDir==-1 & rDir==-1 & rV > bVmax) = legend.id(ismember(legend.name,'reverse_dec')); 
% baseline rev, response rev, response rev within baseline rev range
A(bDir==-1 & rDir==-1 & rV <= bVmax & rV >= bVmin) = legend.id(ismember(legend.name,'no_response')); 
% baseline rev, response pause
A(bDir==-1 & rDir==0) = legend.id(ismember(legend.name,'pause')); 
% baseline rev, response for
A(bDir==-1 & rDir==1) = legend.id(ismember(legend.name,'forward')); 

% baseline pause, response rev
A(bDir==0 & rDir==-1) = legend.id(ismember(legend.name,'reverse')); 
% baseline pause, response pause, response slower than baseline min
A(bDir==0 & rDir==0 & rV < bVmin) = legend.id(ismember(legend.name,'forward_dec')); 
% baseline pause, response pause, response witin baseline range
A(bDir==0 & rDir==0 & rV <= bVmax & rV >= bVmin) = legend.id(ismember(legend.name,'no_response')); 
% baseline pause, response pause, response faster than baseline max
A(bDir==0 & rDir==0 & rV > bVmax) = legend.id(ismember(legend.name,'forward_acc')); 
% baseline pause, response for, response faster than baseline max
A(bDir==0 & rDir==1 & rV > bVmax) = legend.id(ismember(legend.name,'forward_acc')); 
% baseline pause, response for, response witin baseline range
A(bDir==0 & rDir==1 & rV <= bVmax & rV >= bVmin) = legend.id(ismember(legend.name,'no_response')); 
% baseline pause, response for, response slower than baseline max
A(bDir==0 & rDir==1 & rV <bVmin) = legend.id(ismember(legend.name,'forward_dec')); 

% baseline for, response rev
A(bDir==1 & rDir==-1) = legend.id(ismember(legend.name,'reverse')); 
% baseline for, response pause, response faster than baseline max
A(bDir==1 & rDir==0 & rV > bVmax) = legend.id(ismember(legend.name,'forward_acc')); 
% baseline for, response pause, response within baseline range 
A(bDir==1 & rDir==0 & rV <= bVmax & rV >= bVmin) = legend.id(ismember(legend.name,'no_response')); 
% baseline for, response pause, response slower baseline min 
A(bDir==1 & rDir==0 & rV < bVmin) = legend.id(ismember(legend.name,'forward_dec')); 
% baseline for, response for, response faster than baseline max 
A(bDir==1 & rDir==1 & rV > bVmax) = legend.id(ismember(legend.name,'forward_acc')); 
% baseline for, response for, response within baseline range 
A(bDir==1 & rDir==1 & rV <= bVmax & rV >= bVmin) = legend.id(ismember(legend.name,'no_response')); 
% baseline for, response for, response slower than baseline min 
A(bDir==1 & rDir==1 & rV < bVmin) = legend.id(ismember(legend.name,'forward_dec')); 

i = isnan(A);
% find(i)
B = [bDir rDir bVmin bVmax rV A];
if sum(i)~=0
    error('some conditions not accomodated')
end
rDecision = A;


%% convert to text
rDecisionText = legend.name(rDecision);

% %%
% bDirM = repmat(bDir, 1, size(Response, 2));
% 
% 
% 
% %% text
% RT = cell(size(Response)); % declare output array
% RT(dR == 0) = {'baseline'};
% RT(dR == -1 & bDirM == 1 & Response > 0) = {'decelerate forward'};
% RT(dR == -1 & bDirM == 1 & Response < 0) = {'accelerate forward'};
% RT(dR == -1 & bDirM == 1 & Response == 0) = {'pause'};
% RT(dR == -1 & bDirM == -1 & Response < 0) = {'accelerate reverse'};
% RT(dR == -1 & bDirM >= 0 & Response < 0) = {'reversal'};
% RT(dR == 1 & bDirM >= 0 & Response > 0) = {'accelerate forward'};
% RT(dR == 1 & bDirM < 0 & Response < 0) = {'decelerate reverse'};
% RT(dR == 1 & bDirM >= 0 & Response == 0) = {'pause'};
% 
% 
% 
% %% text
% RT = cell(size(Response)); % declare output array
% RT(dR == 0) = {'baseline'};
% RT(dR == -1 & bDirM == 1 & Response > 0) = {'decelerate forward'};
% RT(dR == -1 & bDirM == 1 & Response < 0) = {'accelerate forward'};
% RT(dR == -1 & bDirM == 1 & Response == 0) = {'pause'};
% RT(dR == -1 & bDirM == -1 & Response < 0) = {'accelerate reverse'};
% RT(dR == -1 & bDirM >= 0 & Response < 0) = {'reversal'};
% RT(dR == 1 & bDirM >= 0 & Response > 0) = {'accelerate forward'};
% RT(dR == 1 & bDirM < 0 & Response < 0) = {'decelerate reverse'};
% RT(dR == 1 & bDirM >= 0 & Response == 0) = {'pause'};
% 
% 
% 
% %% numeric
% RTN = nan(size(Response));
% RTN(dR == 0) = 4;
% RTN(dR == -1 & bDirM == 1 & Response > 0) = 5;
% RTN(dR == -1 & bDirM == 1 & Response < 0) = 1;
% RTN(dR == -1 & bDirM == 1 & Response == 0) = 3;
% RTN(dR == -1 & bDirM == -1 & Response < 0) = 2;
% RTN(dR == -1 & bDirM >= 0 & Response < 0) = 7;
% RTN(dR == 1 & bDirM >= 0 & Response > 0) = 1;
% RTN(dR == 1 & bDirM < 0 & Response < 0) = 6;
% RTN(dR == 1 & bDirM >= 0 & Response == 0) = 3;
% % RTN(dR == 1 & baseline_dirM >= 0 & Response == 0) = 3;




%% condense response switch text
% find change data
% a = diff(RTN')';
% 
% b = false(size(RTN,1), size(RTN,2)-1);
% b(a ~=0 & ~isnan(a)) = true;
% 
% c = true(size(RTN,1),1);
% c(isnan(RTN(:,1))) = false;
% 
% b = [c b];
% validata = b;
% 
% % delete repeat response type
% A = RT;
% A(~validata) = {''};
% RTclean = A;
% 
% A = strjoinrows(A,'-');
% % clean 
% A = regexprep(A,'[-]{2,}','-');
% A = regexprep(A,'([-]\>)|(\<[-])','');
% A = regexprep(A,'[-]', ' > ');
% 
% 
% %% save if prompted
% T = tabulate(A);
% 
% T1 = table;
% T1.response = T(:,1);
% T1.n = cell2mat(T(:,2));
% T1.p = cell2mat(T(:,3));




end













































