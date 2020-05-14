function B = DanceM_convert_MWTInd2text(Data,VInd)

B = Data;
fn = B.Properties.VariableNames;
fn1 = fieldnames(VInd);
fn(~ismember(fn,fn1)) = [];
% convert index to variables
for fi = 1:numel(fn)
    B.(fn{fi}) = VInd.(fn{fi})(B.(fn{fi}));
end