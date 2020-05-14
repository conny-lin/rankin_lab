function [pMWTcS,val,pMWTpass] = checkchoroutput(pMWTS,fval,varargin)

displayopt = true;
vararginProcessor;

val = false(size(pMWTS,1),numel(fval));
for mi = 1:numel(pMWTS)
    for fvali = 1:numel(fval)
        if isdir(pMWTS{mi}) == 1
            a = dircontent(pMWTS{mi},fval{fvali});
            if numel(a) > 0; val(mi,fvali) = true; end
        end
    end
end
% reporting
if displayopt
    for fvali = 1:numel(fval)
        fprintf('- %d/%d files do not have %s\n',sum(~val(:,fvali)),size(val,1),fval{fvali});
    end
end

pMWTcS = pMWTS(sum(val,2) ~= numel(fval));
pMWTpass = pMWTS(sum(val,2) == numel(fval));
end