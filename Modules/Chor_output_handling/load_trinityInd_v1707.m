function Tri = load_trinityInd_v1707(pTri,mwtid)

Tri = cell(size(pTri));
Tri_mwtid = Tri;
for i = 1:numel(Tri)
    D = load(pTri{i});
    Tri{i} = table2array(D.Data);
    Tri_mwtid{i} = repmat(mwtid(i),size(D.Data,1),1);
end
Tri = cell2mat(Tri); 
Tri_mwtid = cell2mat(Tri_mwtid);
Tri_mwtid = array2table(Tri_mwtid,'VariableNames',{'mwtid'});
Tri = array2table(Tri,'VariableNames',D.Data.Properties.VariableNames);
Tri = [Tri_mwtid Tri];