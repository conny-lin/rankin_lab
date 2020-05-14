function [pExpf,pzipmwt] = unzipMWT(pIn)
% get paths to folders containing MWT folders
% get all folder content disregard which level it is
[allpaths] = getalldir(pIn); % get paths to all folders under p path

%% get zip files from all dir
pzipf = {};
pExpf = {};
for x = 1:numel(allpaths)
   if isdir(allpaths{x}) ==1;
       [~,p] = dircontent(allpaths{x},'*.zip');
       if isempty(p) ~=1;
           pExpf = [pExpf;allpaths(x)]; % record as it is exp folder
           pzipf = [pzipf;p];
       end
   end
end

%% find zip files that are MWT files
[fn,pzipmwt] = findzipIsMWTf(pzipf); 

%% unzip files
for x = 1:numel(pzipmwt)
    [p,fn] = fileparts(pzipmwt{x,1}); % get home path
    display(sprintf('unzipping files [%s], please wait...',fn));
    cd(p);
    unzip(pzipmwt{x,1},p); % unzip them all
    delete(pzipmwt{x,1}); % delete zip files
    display('done');
end


%% reporting
display(sprintf('%d zipped MWT files found',numel(pzipmwt)));
display(sprintf('in %d Experiment folders contains zipped MWT files',numel(pExpf)));


end
function [zipfnmwt,pzipmwt] = findzipIsMWTf(pzipf)
[a,b] = regexp(pzipf,'\d\d\d\d\d\d\d\d_\d\d\d\d\d\d.zip','match','split');
i = find(not(cellfun(@isempty,a))); % index to MWTfolders
zipfnmwt = {};
pzipmwt = {};
for x = 1:numel(i)
zipfnmwt(x,1) = a{i(x,1),1}(1); % get MWTfoldername
pzipmwt(x,1) = pzipf(i(x,1),1);% get path to MWT folders
end
end



