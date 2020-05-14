function [t] = expsummaryTable(pMWT,varargin)
%% t = expsummaryTable(pMWT,varargin)
% create experiment summary table from pMWT
% VARARGIN:
%     pSave = '/Users/connylin/Dropbox/RL/Dance Output'; % default save output path
%     suffix = ''; % suffix of save file name
%     savetable = false; % default to not save the table in pSave
%     display = true; % default to display table result
%     tableprefix = table; % optional prefix to table

%% DEFAULTS AND PROCESS VARARGIN
pSave = '/Users/connylin/Dropbox/RL/Dance Output'; % default save output path
suffix = ''; % suffix of save file name
savetable = false; % default to not save the table in pSave
display = true; % default to display table result
tableprefix = table; % optional prefix to table
% process varargin
vararginProcessor




%% make table
% parse pMWT into groups
Db = parseMWTinfo(pMWT);
% get exp
enu = unique(Db.expname);
gnu = unique(Db.groupname);
A = nan(numel(enu),numel(gnu));
for ei = 1:numel(enu)
    for gi = 1:numel(gnu)
        en = enu(ei);
        gn = gnu(gi);
        i = ismember(Db.groupname,gn) & ismember(Db.expname,en);
        A(ei,gi) = sum(i);
    end
end
% create table
t = table;
t.expname = enu;
t = [t array2table(A,'VariableNames',gnu)];

% add table prefix
if isempty(tableprefix) == false
    if size(tableprefix,1) ~= numel(enu)
        warning('table prefix not added because row size unmatch');
    else
        t = [tableprefix t];
    end
end

% create sum column
fname = t.Properties.VariableNames;

% get var that's not group name
fn = fname(~ismember(fname,[gnu;{'expname'}]));
a = struct;
for fni = 1:numel(fn)
    fnn = fn{fni};
    d = t.(fnn);
    if isnumeric(d) == 1
        a.(fnn) = NaN;
    elseif iscell(d) == 1
        a.(fnn) = {''};
    end
end

% lable sum under expname
a.expname = {'Total plate by group'};

% sum plate per group
fn = fname(ismember(fname,gnu));
for fni = 1:numel(fn)
    fnn = fn{fni};
    a.(fnn) = sum(t.(fnn));
end
a = struct2table(a);

% join sum row
t = [t;a];

% cal exp plate sum
n1 = size(t,2) - numel(gnu) +1;
n2 = size(t,2);
t.exp_plate_total = sum(table2array(t(:,n1:n2)),2);


%% OPTIONAL
% add suffix (optional)
if isempty(suffix) == 1
    savename = sprintf('%s/exp_summary.csv',pSave);
else
    savename = sprintf('%s/exp_summary %s.csv',pSave,suffix);
end

% save (optional)
if savetable == true
    writetable(t,savename);
end

% display (optional)
if display == true; disp(t); end











