function [pMWTcD,nooutput] = chor_validateoutput(pMWTA,fval)

nooutput = false(size(pMWTA));

% if no output folder, create one
for m = 1:numel(pMWTA)
    if isdir(pMWTA{m}) == 0; 
        mkdir(pMWTA{m}); 
        nooutput(m) = true;
    end
end

display('Checking chor outputs:');
disp(char(fval));

for m = 1:numel(pMWTA)
    for f = 1:numel(fval)
        i = dircontent(pMWTA{m},fval{f});
        if isempty(i) == 1 && nooutput(m) == false
            nooutput(m) = true;
        end
    end
end        
pMWTcD = pMWTA(nooutput);
% pMWTcS = pMWTS(val);
if isempty(pMWTcD) == 1
    disp(sprintf('All %d files contains valid chor output',numel(pMWTA)));
else
    disp('MWT folders below do not contain all chor outputs');
    [~,names] = cellfun(@fileparts,pMWTcD,'UniformOutput',0);
    disp(char(names));
end
    