function [DataGood,DataBad] = validate_trv_N(Data,analysisNLimit)
%% [Data,DataBad] = validate_trv_N(Data,analysisNLimit)
% data is a table containing trv import in column data and mwt path under mwtpath



%% get N
N = nan(size(Data,1),1);
for x = 1:size(Data,1)
    d = Data.data{x};
    N(x) = d.N_alreadyRev(1) + d.N_ForwardOrPause(1) + d.N_Rev(1);
end

ibad = N < analysisNLimit;

%% remove 
DataBad = Data(ibad,:);
DataGood = Data(~ibad,:);

% function [dataval,pMWTval,databad,pMWTbad] = validate_trv_tap(pMWT,Import)
% %% CHECK TAP CONSISTENCY (r20151126)
% 
% %% validate
% if numel(pMWT) ~= numel(Import); error('entries inconsistent'); end
% % create table
% Data = table;
% Data.mwtpath = pMWT;
% Data.d = Import;
% 
% %% get tap size
% [r,c] = cellfun(@size,Data.d,'UniformOutput',0);
% rn = cell2mat(r);
% 
% %% get exp name tap number
% DB = parseMWTinfo(pMWT);
% tapNExpected = DB.tapN;
% 
% %% find trv with inconsistent taps
% i = rn~=tapNExpected;
% pMWTBadTap = Data.mwtpath(i);
% p = pMWTBadTap;
% d = parseMWTinfo(p);
% fprintf('\nPlates with bad taps:\n');
% tabulate(d.groupname)
% 
% %% remove bad taps
% DataBadTap = Data(i,:);
% Data(i,:) = [];
% 
% %% report
% p = Data.mwtpath;
% d = parseMWTinfo(p);
% fprintf('\nPlates with validated taps:\n');
% tabulate(d.groupname)
% 
% %% output
% pMWTval = Data.mwtpath;
% dataval = Data.d;
% pMWTbad = DataBadTap.mwtpath;
% databad = DataBadTap.d;



% %% exclude data with nan frequency data
% i(any(isnan(MWTSet.Data.Raw.Y.RevFreq))) = false;
% % reporting
% mwtfn = MWTSet.Data.Raw.MWTfn(~i);
% disp(char(mwtfn));
% fprintf('\n\nRemove from analysis: <%d Total worms at first tap + nan value for freq\n',analysisNLimit);
% D = MWTSet.Data.Raw;
% B = struct;
% Bad = struct;
% type = fieldnames(D);
% for ti = 1:numel(type)
%     d = D.(type{ti});
%     if isstruct(d) == 1
%         fnames = fieldnames(d);
%         for x =1:numel(fnames)
%             B.(type{ti}).(fnames{x}) = d.(fnames{x})(:,i);
%             Bad.(type{ti}).(fnames{x}) = d.(fnames{x})(:,~i);
%         end
%     elseif size(d,1) == numel(i)
%         B.(type{ti}) = d(i,:);
%         Bad.(type{ti}) = d(~i,:);
%     elseif size(d,2) == numel(i)
%         B.(type{ti}) = d(:,i);
%         Bad.(type{ti}) = d(:,~i);
%     else
%         error('type not accomodated')
%     end
% end
% MWTSet.Data.Raw = B;
% MWTSet.Data.Raw_Excluded = Bad;
% % report
% p = B.pMWT;
% t = parseMWTinfo(p);
% fg = t.groupname;
% if numel(fg) == 0
%     fprintf('\nNo plates qualified, abort\n');
%     cd(pSave);
%     fid = fopen('No plates qualified for this analysis.txt','w');
%     fclose(fid);
%     return
% else
%     fprintf('\nValide Plates:\n');
%     tabulate(fg)
% end