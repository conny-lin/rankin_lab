function T = translate_rx2alcohol(T)

% -------------------------------------------------------------------------
%                   INPUT VALIDATION | 20170822
% -------------------------------------------------------------------------
if ~istable(T) % check if table
    error('function not accomodating non tables')
else
    % check if has field names
    if ~ismember('rx',fieldnames(T))
        if ~ismember('dose',fieldnames(T))
            error('no fieldname rx or dose')
        end
    end
end

% -------------------------------------------------------------------------
%                   PROCESS | 20170822
% -------------------------------------------------------------------------
%%
if ismember('rx',fieldnames(T))
  i = ismember(T.rx ,'NA');
  T.rx(i) = {'0mM'};
  T.dose = T.rx;
end
if ~ismember('dose',fieldnames(T))
   T.dose = T.rx; 
end
