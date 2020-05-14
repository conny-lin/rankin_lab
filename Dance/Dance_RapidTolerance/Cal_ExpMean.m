

%% mean per experiment
D = MWTSet.Data_Plate;
D = innerjoin(MWTSet.Info.MWTDB,D,'Keys','mwtid');
% create exp-groupname combo
a = [D.expname D.groupname_short];
b = cell(size(a,1),1);
for x = 1:size(a,1)
    b{x} = strjoin(a(x,:),' x ');
end
D.expXgname = b;


% find means by exp
grp = D.expXgname;
gnameu = unique(grp); 
% make header
grpheader = table;  
grpheader.expXgname = gnameu;
a = regexpcellout(gnameu,' x ','split');
grpheader.expname = a(:,1);
a = a(:,2);
a = regexpcellout(a,' ','split');
grpheader.strain = a(:,1);
grpheader.predose = a(:,2);
grpheader.postdose = a(:,3);
% declare struct array
A = struct;
for msri = 1:numel(msrlist)
    T = grpstatsTable(D.(msrlist{msri}),grp,'grpheader',grpheader,'gnameu',gnameu);
    [i,j] = ismember(T.gnameu,gnameu);
    T = [grpheader(j,:) T];
    T.gnameu = [];
    
    A.(msrlist{msri}) = T;
end
MWTSet.Data_Exp = A;

%% export mean by experiment
% D = ;
for msri = 1:numel(msrlist)
    D = MWTSet.Data_Exp.(msrlist{msri});
    grp = strjoinrows(D(:,{'strain','predose','postdose'}));
    T = grpstatsTable(D.mean,grp);
    a = regexpcellout(T.gnameu,' ','split');
    
    A = table;
    A.strains = a(:,1);
    A.predose = a(:,2);
    A.postdose = a(:,3);
    T = [A T];
    T.gnameu = [];
    MWTSet.Data_Group_byExp.(msrlist{msri}) = T;
    filename = sprintf('%s/%s Nexp.csv',pSave,msrlist{msri});
    filename = sprintf('%s/%s Nexp.csv',pSave,msrlist{msri});

    writetable(T,filename)
end

%% export as graphic output

for msri= 1:numel(msrlist)
    D = MWTSet.Data_Group_byExp.(msrlist{msri});

    prefix = uniqueCellrows(D(:,{'predose','postdose'}));
    cond = strjoinrows(D(:,{'predose','postdose'}));
    condu = unique(cond);
    % make col names
    su = unique(D.strains);
    su = [su(ismember(su,'N2')); su(~ismember(su,'N2'))];
    var = [su strjoinrows([su repmat({'se'},numel(su),1)],'_')];
    %% assign
    A = nan(numel(condu),numel(var));
    for ci = 1:numel(condu)
        i = ismember(cond,condu{ci});
        s = D.strains(i);
        m = D.mean(i);
        se = D.se(i);
        [i,j] = ismember(su, s);
        A(ci,:) = [m(j(i))' se(j(i))'];
    end
    T = [prefix array2table(A,'VariableNames',var)];
    filename = sprintf('%s/%s Nexp graphsetup.csv',pSave,msrlist{msri});
    writetable(T,filename);
end






