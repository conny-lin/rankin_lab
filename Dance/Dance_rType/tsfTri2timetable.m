function S = tsfTri2timetable(Tri,varnames,timename,timeint)


timestep = 1/timeint;

Tri.(timename) = seconds(Tri.(timename));
S = struct;
for vi = 1:numel(varnames)
    v = varnames{vi};
    D = Tri(:,{'ids',timename,v});
    S.(v).timetable_raw = D;
    
    D.(timename) = D.(timename)*timestep;
    B = table2timetable(D,'RowTimes',timename);
    B = sortrows(B,{timename});
    B = unstack(B,v,'ids','AggregationFunction',@nanmean);
    B = retime(B,'secondly',@nanmean);
    B.(timename) = B.(timename)./timestep;
    S.(v).timetable = B;
    S.(v).array = table2array(B);
end

S = remove_nan_data(S,'timeparts','timename',timename);
