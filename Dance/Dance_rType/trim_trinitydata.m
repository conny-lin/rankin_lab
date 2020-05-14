function Tri = trim_trinitydata(Tri,starttime,finishtime,varnames,varargin)

removenan = false;
% varnames = {'id','ids','mwtid','wormid','time','speed','bias','tap'};
% delete data outside of assay times
Tri(Tri.time < starttime | Tri.time > finishtime,:) = [];
Tri.id = Tri.wormid./10000 + Tri.mwtid; % calculate id
a = strjoinrows([cellfunexpr(Tri.mwtid,'w') num2cellstr(Tri.mwtid) num2cellstr(Tri.wormid)],'_'); % create id string
Tri.ids = a;
Tri = Tri(:,varnames);
%--------------------------------------------------------------
if removenan
    Tri(isnan(Tri.bias) | isnan(Tri.speed),:) = [];
end

