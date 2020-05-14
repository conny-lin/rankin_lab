function [MWTSet,textout,DS] = Stats_HabCurveRMANOVA_v1707(MWTSet,pSave,varargin)
% -------------------------------------------------------------------------
%% default
% -------------------------------------------------------------------------
% alpha = 0.05;
% pvlimit = 0.001;
prefix = '';
suffix = '';
msrlist = {'RevFreq','RevSpeed','RevDur'};
%--------------------------------------------------------------------------
%% varargin 
% -------------------------------------------------------------------------
vararginProcessor
%--------------------------------------------------------------------------
%% prepare data 
% -------------------------------------------------------------------------
M = MWTSet.MWTDB;
Data = MWTSet.Raw;
A = innerjoin(Data, M(:,{'mwtid','groupname','strain','rx'}));
factors = {};
if numel(unique(A.strain))>1
    factors = {'strain'};
end
if any(regexpcellout(A.rx,'mM'))
    A.rx(regexpcellout(A.rx,'NA')) = {'0mM'};
    A.dose = regexpcellout(A.rx,'\d{1,}(?=mM)','match');
    factors = [factors {'dose'}];
end
if isempty(factors)
    factors = {'groupname'};
end
Data = A;
%--------------------------------------------------------------------------
%% anova | 20170816
% -------------------------------------------------------------------------
filesavename = fullfile(pSave,sprintf('%sRMANOVA%s.txt',prefix,suffix));
fid = fopen(filesavename,'w');
for msri = 1:numel(msrlist)
    msr = msrlist{msri};
    % anova
    [textout,DS] = anovarm_std(Data,'tap',factors,'mwtid',msr);
    % write
    if msri==1
        fprintf(fid,'----- %s -----\n%s\n',msr,textout);
    else
        fprintf(fid,'\n----- %s -----\n%s\n',msr,textout);
    end
    % store in data
    MWTSet.Stats.(msr).Curve = DS;  
end
fclose(fid);
% -------------------------------------------------------------------------
%% anova by strain | 20170816
% -------------------------------------------------------------------------
for msri = 1:numel(msrlist)
    msr = msrlist{msri}; % get measure name
    D = MWTSet.Stats.(msr).Curve.Data;    % get data
    strainu =unique(D.strain);    % get strain names
    for si = 1:numel(strainu)
        DStrain = D(ismember(D.strain,strainu{si}),:);
        [textout,ANOVAS] = anovarm_std(DStrain,'tap',{'dose'},'mwtid',msr); % anova
        MWTSet.Stats.(msr).Curve_by_strain = ANOVAS; % store in MWTSet
    end
end
% -------------------------------------------------------------------------
















