function D = convert_import2table(Data,legend,varargin)

%% varargin input varables default
msr = legend;
vararginProcessor


%% validate legend
if size(Data,2) ~= numel(legend); error('incorrect legend'); end
if iscell(Data)
    D=cell2mat(Data);
else
    D=Data;
end

[i,j] = ismember(msr,legend);
D = D(:,j(i));
D = array2table(D,'VariableNames',msr(i));
