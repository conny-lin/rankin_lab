function str = gen_Nstring(n,varargin)
%% generate n string


%% DEFAULTS & VARARGIN
nInput = 1;
sampleName = 'plate';

% varargin processer
if nargin > nInput
    A = varargin;
    if isinteger(numel(A)/2) == 1
        error('variable entries incorrect')
    end
    callName = A(1:2:numel(A));
    % check if call names are strings
    for x = 1:numel(callName)
       if ischar(callName{x}) == 0
           error('varargin incorrect');
       end
    end
    % get values
    setValue = A(2:2:numel(A));
    for x = 1:numel(callName)
        % determine input type
        if eval(sprintf('ischar(%s)',callName{x})) == 1
            eval(sprintf('%s = ''%s'';',callName{x},setValue{x}));
        elseif eval(sprintf('isnumeric(%s)',callName{x}))
            eval(sprintf('%s = %d;',callName{x},cell2mat(setValue(x))));
        elseif eval(sprintf('iscell(%s)',callName{x}))
            a = setValue{x};
            eval(sprintf('%s = a;',callName{x}));
        else
            error('need coding');
        end
    end
end



%% produce strin
str = sprintf('N(%s) = ',sampleName);
for x = 1:numel(n)
    if x ~=numel(n);
        str = sprintf('%s%d, ',str,n(x));
    else
        str = sprintf('%s%d',str,n(x));
    end
end



end