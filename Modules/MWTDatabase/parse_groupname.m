function varargout = parse_groupname(gn,query,varargin)


if nargin == 1
    query = {'strain','rx'};
end
vararginProcessor;

%%
if size(gn,1)==1
    gn = gn';
elseif size(gn,2)==1
%     gn = gn';
end

%%
varargout = cell(size(query));
for x = 1:numel(query)
    qname = query{x};
    switch qname
        case 'strain'
            varargout{x} = regexpcellout(gn,'\<[A-Z]{1,}\d{1,}','match');
        case 'rx'
            a = regexpcellout(gn,'_','split');
            varargout{x} = regexpcellout(gn,'(?<=\<[A-Z]{1,}\d{1,}_)\w{1,}','match');
     
    
    end
end