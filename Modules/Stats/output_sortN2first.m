function [T1,ind] = output_sortN2first(T,varargin)

tblcolname = 'gname';
vararginProcessor


if istable(T)
    a = find(regexpcellout(T.(tblcolname),'N2'));
    b = find(~regexpcellout(T.(tblcolname),'N2'));
    ind = [a;b];
    T1 = T(ind,:);
elseif size(T,1)==1 || size(T,2)==1
    a = find(regexpcellout(T,'N2'));
    b = find(~regexpcellout(T,'N2'));
    ind = [a;b];
    T1 = T(ind);
else
    error('code for this scenario')
end
