msrlist = {'speed','curve'};
S = innerjoin(MWTDB,MWTSet.Data_Plate,'Keys','mwtid');

strainu = unique(S.strain);
for msri = 1:numel(msrlist)
    % get measure name
    msr = msrlist{msri};
    % current data
    Y = S.(msr);
    % multi-factorial anova
    gvarnames = {'strain','predose','postdose'};
    G = {S.strain, S.predose, S.postdose};
    anovan_std(Y,G,gvarnames,pSave,'export','textsave','suffix',msrlist{msri});

    %% within strain effect
    for si = 1:numel(strainu)
        strainame = strainu{si};
        i = ismember(S.strain,strainame);
        Y1 = Y(i);
        G1 = {S.predose(i), S.postdose(i)};
        gvarnames = {'predose','postdose'};
        anovan_std(Y1,G1,gvarnames, pSave,'prefix', msrlist{msri},'suffix',strainu{si},'export','text');
    end
end