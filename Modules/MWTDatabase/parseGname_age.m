
function B = parseGname_age(G)

%% input can be groupname or MWTDB
%%
if istable(G)
   gn = G.groupname;
elseif iscell(G) && width(G)==1
   gn = G;
end

%%
B = cell(size(gn));
b = regexpcellout(gn,'(?<=[_])\d{1,}(?=d)','match');
i = ~cellfun(@isempty,b);
B(i) = b(i);
B(~i) = {'4'};


if istable(G)
   G.age_test = B;
      B = G;

end