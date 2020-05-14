function [A] = celltakeout(B,varargin)

%% if no option selection, examine input B
switch nargin
    case 1
        
       %% 2 level cell
        if iscell(B{1}) ==1
            i = cell2mat(cellfun(@size,B,'UniformOutput',0));
            onecolume = sum(i(:,2)==1) == numel(B); % see if only one colume
            multirow = sum(i(:,1)>1) ~= 0; % see if multiple rows
            if onecolume == 1 %onecolume == 1 && multirow == 1
                option = 'multirow';
            elseif onecolume == 0 && multirow == 0
                option = 'split';
            end
        end
        
        %% 1 level cell
        if iscell(B{1}) == 0 && iscell(B) ==1
            if size(B,2) >1
                % see if only one colume
                onecolume = 0;
            else 
                onecolume = 1;
            end
            if size(B,1) > 1
                % see if multiple rows
                multirow = 1;
            else 
                multirow = 0;

            end

            if onecolume == 1 && multirow == 1
                option = 'multirow';
            end
            
            if onecolume == 0 && multirow == 0
                option = 'split';
            end
        end
        
    case 2
        option = varargin{1};
       
end

 %% if option is given
switch option
    case 'onecolmultirow'

    case 'split'
        A = {};
        for x = 1:numel(B); 
            col = size(B{x},2); 
            A(x,1:col) = B{x}; 
        end
    case 'multirow'
        A = {}; 
        for x = 1:numel(B); 
            A = [A;B{x}]; 
        end
    case 'singlerow'
        A = {};
        for x = 1:numel(B);
            if isempty(B{x})==0; A(x,1) = B{x,1};
            else A(x,1) = {''}; end
        end
    case'match'
        A = {};
        for x = 1:numel(B);
            if isempty(B{x})==0; A(x,1) = B{x,1};
            else A(x,1) = {''}; end
        end
    case 'singlenumber'
        A = [];
        for x = 1:numel(B)
            if isempty(B{x})==0; A(x,1) = B{x,1};
            else A(x,1) = 0; end
        end
    case 'logical'
         A = [];
        for x = 1:numel(B)
            if isempty(B{x})==0; A(x,1) = B{x,1};
            else A(x,1) = 0; end
        end
        A = logical(A);

    otherwise
        error('no %s function',option);
end

end






