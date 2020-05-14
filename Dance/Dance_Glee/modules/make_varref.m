function [VInd,D,varref,GroupVarName] = make_varref(D,VarName,VInd)

%%
vn = D.Properties.VariableNames;

%% create VarName combintation
% check if input has VarName variables
if sum(ismember(VarName,vn)) ~= numel(VarName)
    error('some VarName input are not found');

else
    % make VarName table
    t = table;
    for x = 1:numel(VarName)
       t.(VarName{x}) = D.(VarName{x});
    end
    varref = t;
    varcombo = unique(t,'rows');
end


%% make combined group name and create index if VarName > 1

if numel(VarName) > 1
    % find actual names from VInd
    vn = varcombo.Properties.VariableNames;
    A = varcombo;
    for x = 1:numel(vn)
        v = vn{x};
        A.(v) = VInd.(v)(A.(v));
    end
    % join group name
    B = cell(size(A,1),1);
    A = table2array(A);
    for x = 1:size(A,1)
        B{x} = strjoin(A(x,:),'*');
    end
    % make new combined group name
    for x = 1:numel(VarName)
        if x == 1
            a = VarName{x};
        else
            b = VarName{x};
            a = sprintf('%s_%s',a,b);
        end
    end
    VInd.(a) = B;
    
    %% make combined varref and varcombo
    [~,varref] = ismember(varref,varcombo,'rows');
    D.(a) = varref;
    GroupVarName = a;
    
    
end  
    