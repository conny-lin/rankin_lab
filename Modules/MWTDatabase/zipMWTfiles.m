function zipMWTfiles(pFolder,type) 
%% zip MWT files
% pFolder = '/Volumes/FLAME/Conny_Lin/20150910C_CL_100s30x10s10s_slo1rescue1';

%% process input
if nargin < 2
    type = 'expfolder';
end
%% get whats in the mwt folder
switch type
    case 'expfolder'
        
        [~,~,~,pf] = dircontent(pFolder);
        for x = 1:numel(pf)
            [~,~,mwtf,pmwt] = dircontent(pf{x});
            for y = 1:numel(mwtf)
                zip(pmwt{y},pmwt{y});
                rmdir(pmwt{y},'s');
            end
        end
        
    case 'MWTfolder'
        [~,~,mwtf,pmwt] = dircontent(pf{x});
        for y = 1:numel(mwtf)
            zip(pmwt{y},pmwt{y});
            rmdir(pmwt{y},'s');
        end
        
        
end
display 'DONE';
end

%% SUBFUN

function [varargout] = dircontent(p,varargin)
% a = full dir, can be folder or files, b = path to all files, 
% c = only folders, d = paths to only folders

%% get dir
switch nargin
    case 0
        varargout{1} = [];
        return
    
    case 1
        % get dir
        [a,b,c,d] = dirc(p);
        varargout{1} = a;
        varargout{2} = b;
        varargout{3} = c;
        varargout{4} = d;
        return
    
    
    case 2
        ext = varargin{1};
        % dircontnetext
        cd(p); % go to directory
        a = {};
        a = dir(ext); % list content
        a = {a.name}'; % extract folder names only
        a(ismember(a,{'.','..','.DS_Store'})) = []; 
        b = {};
        for x = 1:size(a,1); % for all files 
            b{x,1} = [p,'/',a{x,1}]; % make path for files
        end
         varargout{1} = a;
         varargout{2} = b;
         return
    
         
    case 3
        if strcmp(varargin{1},'Option')
            option = varargin{2};
        end
        
end


switch option
    case 'MWT'
        p1 = regexp(genpath(p),':','split')';
        [~,fn2] = cellfun(@fileparts,p1,'UniformOutput',0);
        i = ismember(fn2,{''});
        fn2(i) = []; 
        p1(i) = [];
        i = regexpcellout(fn2,'\<\d{8}[_]\d{6}\>');
        varargout{1} = fn2(i);
        varargout{2} = p1(i);
    case 'MWTall'
        p1 = regexp(genpath(p),':','split')';
        [~,fn2] = cellfun(@fileparts,p1,'UniformOutput',0);
        i = ismember(fn2,{''});
        fn2(i) = []; 
        p1(i) = [];
        i = regexpcellout(fn2,'\<\d{8}[_]\d{6}');
        varargout{1} = fn2(i);
        varargout{2} = p1(i);

    otherwise
        error 'No such option';
end
end




%% Subfunction
function [a,b,c,d] = dirc(p)
% a = full dir, can be folder or files, b = path to all files, 
% c = only folders, d = paths to only folders
cd(p); % go to directory
a = {}; % create cell array for output
a = dir; % list content
a = {a.name}'; % extract folder names only
a(ismember(a,{'.','..','.DS_Store'})) = []; 
b = {};
c = {};
d = {};
for x = 1:size(a,1); % for all files 
    b(x,1) = {strcat(p,'/',a{x,1})}; % make path for files
    if isdir(b{x,1}) ==1; % if a path is a folder
        c(end+1,1) = a(x,1); % record name in cell array b
        d(end+1,1) = b(x,1); % create paths for folders
    else
    end
end
end