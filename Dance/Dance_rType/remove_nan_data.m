function S = remove_nan_data(S,removenan,varargin)

timename = 'Time';
vararginProcessor;
%% remove nan
i = false;
varnames = fieldnames(S);
switch removenan
    case 'any'
        i = any(isnan(S.bias.array) & isnan(S.speed.array));

    case 'all'
        i = sum(isnan(S.bias.array)) == size(S.bias.array,1);
    case 'timeparts'
        t1 = S.bias.timetable.(timename)<0;
        t2 = S.bias.timetable.(timename)>0;
        B = S.bias.array(t1,:);
        R = S.bias.array(t2,:);
        i = sum(isnan(B)) == size(B,1) | sum(isnan(R)) == size(R,1);
        
end

%%

% if any(i)
    m = {'array','timetable'};
    for mi = 1:numel(m)
        for vi = 1:numel(varnames)
            S.(varnames{vi}).(m{mi})(:,i) = [];
        end
    end
% end
