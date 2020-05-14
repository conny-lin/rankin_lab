[Data,~,Legend2] = import_drunkposture2_dat(pMWT,'array');
if numel(pMWT) ~=numel(Data); error('import number incorrect'); end
Data = array2table(cell2mat(Data),'VariableNames',Legend2);
% remove data outside of assay time
Data(Data.time < timeStartSet | Data.time > timeEndSet,:) = [];
% store in MWTSet
MWTSet.Import.drunkposture2 = Data;
MWTSet.Info.MWTDB = MWTDB;
MWTSet.Info.legend.drunkposture2 = Legend2;
% clear memory
clear gnn s a b c strain expr e i
