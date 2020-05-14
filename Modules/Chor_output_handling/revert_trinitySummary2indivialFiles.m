function revert_trinitySummary2indivialFiles(pMWTp,overwrite,varargin)

displayopt = false;
if nargin < 2
    overwrite = false;
end
vararginProcessor;
%% see if exist
prefix = 'trinity_wormid';
fn = dircontent(pMWTp,[prefix,'*']);
wormid_exist = [];
if ~isempty(fn)
    searchstr = sprintf('(?<=%s)\\d{1,}(?=[_])',prefix);
    wormid_exist = cellfun(@str2num,regexpcellout(fn,searchstr,'match'));
end


%% get data
Trinity = load_Trinity_v1707(pMWTp);
D = extract_trinity2table(Trinity); % convert to table
A = grpstats(D,'wormid',{'min','max','numel'},'DataVars',{'time'}); % extract times by worms

if ~overwrite
   i = ismember(A.wormid, wormid_exist);
   if sum(i) == numel(i); 'break'; return; end
   A(i,:) = []; % delete worms already done
end

%% write data
if displayopt; fprintf('*** Converting Trinity summary to inidividual files ***\n');
for wormi = 1:size(A,1)
    if displayopt; loopreporter(wormi,'',20,size(A,1)); end
    wormid = A.wormid(wormi);
    t1 = A.min_time(wormi);
    t2 = A.max_time(wormi);
    n = A.numel_time(wormi);
    savepath = fullfile(pMWTp,sprintf('%s%d_%.0fs_%.0fs_N%d.mat',prefix,wormid,t1,t2,n));
    Data = D(D.wormid==wormid,:); % get data
    save(savepath,'Data');
end

end