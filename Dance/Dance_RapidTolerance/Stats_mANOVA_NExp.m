msrlist = {'speed','curve'};
Data = MWTSet.Data_Exp;
for msri = 1:numel(msrlist)
    % get measure name
    msr = msrlist{msri};
    D = Data.(msr);
    strainu = unique(D.strain);

    % current data
    Y = D.mean;
    % multi-factorial anova
    gvarnames = {'strain','predose','postdose'};
    G = {D.strain, D.predose, D.postdose};
    prefix = sprintf('%s NExp', msrlist{msri});
    anovan_std(Y,G,gvarnames,pSave,'prefix',prefix,'export','text');
    
    %% within strain effect
    for si = 1:numel(strainu)
        strainame = strainu{si};
        i = ismember(D.strain,strainame);
        Y1 = Y(i);
        G1 = {D.predose(i), D.postdose(i)};
        gvarnames = {'predose','postdose'};
        suffix = sprintf('%s NExp',strainu{si});
        anovan_std(Y1,G1,gvarnames, pSave, 'prefix',msrlist{msri},'suffix',suffix,'export','text');
    end
end