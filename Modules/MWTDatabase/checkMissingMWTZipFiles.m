function [pMWTnew,pMWTdup] = checkMissingMWTZipFiles(pCheck,MWTdb)


%% find available files
pMWTzip = getMWTzipped(pCheck);
MWTDBzip = parseMWTinfo(pMWTzip);

%% compare
mwtf1 = MWTdb.mwtname;
mwtf2 = MWTDBzip.mwtname;
i = ~ismember(mwtf1,mwtf2);
mwtfm = mwtf1(i);
disp(mwtfm);
pMWTnew = MWTdb.mwtpath(i);
disp(unique(MWTdb.expname(i)))
%% get pmwtpath
MWTndb = parseMWTinfo(pMWTnew);

%% compare if missing paths has exp name folder and group folder matched
e1 = unique(MWTndb.expname);
e2 = unique(MWTDBzip.expname);
e1m = e1(ismember(e1,e2));
display 'these expname already exist';
disp(e1m)
%% get paths with existing expname
pMWTdup = MWTndb.mwtpath(ismember(MWTndb.expname,e1m));




%% diplay
for x = 1:numel(e1m)
    % get mwt name
    [~,mn] = cellfun(@fileparts,pMWTnew(regexpcellout(pMWTnew,e1m{x})),'UniformOutput',0);
    % display
    fprintf('\n%s\n',e1m{x});
    disp(char(mn));
end

display 'these expname do not exist';
disp(e1(~ismember(e1,e2)))
return
%% find closer match for exp date
% e1md = regexpcellout(e1m,'\<\d{8}','match');
%%
% e = MWTDBzip.expname(ismember(MWTDBzip.exp_date,cellfun(@str2num,e1md)));
% unique(e)
