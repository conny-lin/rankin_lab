function [MWTSet] = DanceM_exportcsv(MWTSet,varargin)


%% default
% displayopt = false;
% overwrite = false;
dancename = MWTSet.ProgramName;

%% varargin
vararginProcessor;


%% export as graphic output
switch dancename
    case 'Dance_RapidTolerance_v1707'
     otherwise
        error('no chor type determined for %s',dancename);   
end


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