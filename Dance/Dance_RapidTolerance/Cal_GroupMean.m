msrlist = {'speed','curve'};
%%
A = struct;
grp = S.groupname_short;
gnameu = unique(grp);
% make header
grpheader = table;  
a = regexpcellout(gnameu,' ','split');
grpheader.strains = a(:,1);
grpheader.predose = a(:,2);
grpheader.postdose = a(:,3);
for msri= 1:numel(msrlist)
    T = grpstatsTable(S.(msrlist{msri}), grp,'grpheader',grpheader,'gnameu',gnameu);
    %%
    a = regexpcellout(T.gnameu,' ','split');
    T1 = table;
    T1.strains = a(:,1);
    T1.predose = a(:,2);
    T1.postdose = a(:,3);
    T = [T T1];
    T.gnameu = [];
    
    %%
    filename = sprintf('%s/%s Nplate.csv',pSave,msrlist{msri});
    writetable(T,filename)
    A.(msrlist{msri}) = T;
end
MWTSet.Data_Group = A;

%%

save(fullfile(pSave,'test.mat'),'A');

%% export as graphic output
for msri= 1:numel(msrlist)
    D = MWTSet.Data_Group.(msrlist{msri});
    prefix = uniqueCellrows(D(:,{'predose','postdose'}));
    cond = strjoinrows(D(:,{'predose','postdose'}));
    condu = unique(cond);
    % make col names
    su = unique(D.strains);
    su = [su(ismember(su,'N2')); su(~ismember(su,'N2'))];
    a = strjoinrows([su repmat({'se'},numel(su),1)],'_');
    var = [su; a];
    %% assign
    A = nan(numel(condu),numel(var));
    for ci = 1:numel(condu)
        i = ismember(cond,condu{ci});
        s = D.strains(i);
        m = D.mean(i);
        se = D.se(i);
        [i,j] = ismember(su, s);
        a = [m(j(i))' se(j(i))'];
        A(ci,:) = a;
    end
    T = [prefix array2table(A,'VariableNames',var)];
    filename = sprintf('%s/%s Nplate graphsetup.csv',pSave,msrlist{msri});
    writetable(T,filename);
end











