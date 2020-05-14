function arraySorted = sortN2first(groupnamelist, array)
    

i = regexpcellout(groupnamelist,'N2');
i = [find(i) ;find(~i)];
% check array type
[rn,cn] = size(array);
if rn == numel(groupnamelist)
    arraySorted = array(i,:);
elseif cn == numel(groupnamelist)
    arraySorted = array(:,i);
else
    error('array input size does not match group number');
end
