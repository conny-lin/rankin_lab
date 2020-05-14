function T = saveGraphstruct2table(GData,pSave)
%%
% GData = G.AccProb;


%% transform G struct to graph table


c1 = GData.N(:,1);
c2 = GData.X(:,1);
c3_6 = GData.Y;
c7_10 = GData.E;
GD = [c1 c2 c3_6 c7_10];

% make column names
gnames = GData.groupname;
leg = {'Y','E'};
legtitle = {'N';'Time'};

for il = 1:numel(leg)
    a = cellfunexpr(gnames,leg{il});
    a = [gnames a];
    b = cell(size(gnames));
    for i = 1:size(a,1)
        b(i) = strjoinrows(a(i,:),'_');
    end
    legtitle = [legtitle;b];
end

%%
T = array2table(GD,'VariableNames',legtitle);

writetable(T,pSave);