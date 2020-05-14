function [MWTSet,Data] = DanceM_datatransform(MWTSet,varargin)

%% default
dancename = MWTSet.ProgramName;

%% varargin
vararginProcessor


%% determine variables wanted from dance name
switch dancename
    case 'Dance_RapidTolerance_v1707'
        v_wanted = {'mwtid','time','goodnumber','speed','bias','curve'};
        transform_wanted = {'group_by_exp'};
    otherwise
        error('no chor type determined for %s',dancename);

end

%% check if all variable exists
D = MWTSet.Raw;
va = D.Properties.VariableNames;
vi = intersect(va,v_wanted);
if numel(vi)~=numel(v_wanted)
    error('raw data does not include certain variables needed');
end
% trim raw data
D = D(:,v_wanted);


%% process
MWTDB = MWTSet.MWTDB;
Data = struct;
for twi = 1:numel(transform_wanted)
    twname = transform_wanted{twi};
    switch twname
        case 'group_by_exp'
            M = MWTDB(:,{'mwtid','expname','groupname_short'});
            M.groupname = M.groupname_short;
            M.G = strjoinrows([M.expname, M.groupname],'@');
            D = innerjoin(M,D);
            A = grpstats(D,{'expname','groupname'},{'numel','mean','sem'},'DataVars',{'speed','curve'});
            Data.(twname) = A;
            return
    end
end
MWTSet.Data = Data;