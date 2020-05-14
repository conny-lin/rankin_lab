function export_drunkposture2_bygroup(D,pSaveA,varargin)

%% defaults
msru = fieldnames(D);
gname = 'groupname';
prefix = '';
vararginProcessor
if isdir(pSaveA) == 0; mkdir(pSaveA); end

%% transform

for mi = 1:numel(msru)
    msr = msru{mi};
%%     fnu = fieldnames(D.(msr));
    d = D.(msr);
    T = table;
    T.(gname) = d.(gname);
    
    t = d.mean';
    t = array2table(t);
    T = [T t];
    SE = d.SE';
    t = array2table(SE);
    T = [T t];
    
    %%
    writetable(T,sprintf('%s/%s%s.csv',pSaveA,prefix,msr));
end