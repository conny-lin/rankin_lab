function T = split_genotype_info(genotypes)
% genotypes = S.genotype;
a = regexpcellout(genotypes,'[(]|[)]|[;]|[[]|[]]','split');

%% go through each strain
varnames = {'genotype','promoter','gene','allele','chromosome','array','others'};
A = cell(size(genotypes,1),numel(varnames));
T = cell2table(A,'VariableNames',varnames);
startn = 1;

for ai = startn:size(a,1)
    % genotype
    genotype = genotypes{ai};
    T.genotype(ai) = {genotype};
    % get info
    b = a(ai,:)';
    b(cellfun(@isempty,b)) = [];
    for vi = 2:numel(varnames)
        if ~isempty(b)
            name = varnames{vi};
            switch name
                case 'gene'
                    searchterm = '(^[a-zA-Z]{1,}-\d{1,}$)|(^wildtype$)|(^[A-Z]\d+[A-Z]+\d+[.]\d+)';
                case 'chromosome'
                    searchterm = '^((I){2,3})|^I$(^IV$)|^I$|X|V$';
                case 'promoter'
                    searchterm = '^(P|p)[a-zA-Z]{1,}-\d{1,}[:]{1,}[a-zA-Z]{1,}-\d{1,}';
                case 'array'
                    searchterm = '(^(EX)|(Ex)([a-z]{1,}|[0-9]{1,}))|([a-z]{1,}Is[0-9]{1,})';
                case 'allele'
                    searchterm = '^[a-z]{1,}\d{1,}|(hawaiian)|(bristol)';
            end
            switch name 
                case 'others'
                    c = strjoinrows(b',', ');
                    
                otherwise
                    % process
                    c = regexpcellout(b,searchterm,'match');
                    i = regexpcellout(b,searchterm);
                    b(i) = [];

                    c(cellfun(@isempty,c)) = [];
                    if numel(c)>1
                        c = strjoinrows(c',', ');
                    end

            end
            if ~isempty(c)
                T.(name)(ai) = c;
            end
        end

    end
end


%% go through each strain
% varnames = {'genotype','gene','promoter','allele','chromosome','array'};
% A = cell(size(genotypes,1),numel(varnames));
% T = cell2table(A,'VariableNames',varnames);
% startn = 1;
% genotype_parts = regexpcellout(genotypes,'[;]','split');
% 
% 
% for ai = startn:size(genotype_parts,1)
%     % genotype
%     genotype = genotypes{ai};
%     T.genotype(ai) = {genotype};
%     % get info
%     b = genotype_parts(ai,:)';
%     b(cellfun(@isempty,b)) = [];
%     for vi = 2:numel(varnames)
%         if ~isempty(b)
%             name = varnames{vi};
%             switch name
%                 case 'gene'
%                     searchterm = '(^[a-zA-Z]+(-)\d+(?=[(])|(^wildtype$)|(^[A-Z]\d+[A-Z]+\d+[.]\d+)';
%                 case 'chromosome'
%                     searchterm = '^((I){2,3})|^I$(^IV$)|^I$|X|V$';
%                 case 'promoter'
%                     searchterm = '^(P|p)[a-zA-Z]{1,}-\d{1,}[:]{1,}[a-zA-Z]{1,}-\d{1,}';
%                 case 'array'
%                     searchterm = '(^(EX)|(Ex)([a-z]{1,}|[0-9]{1,}))|([a-z]{1,}Is[0-9]{1,})';
%                 case 'allele'
%                     searchterm = '^[a-z]{1,}\d{1,}|(hawaiian)|(bristol)';
%             end
% 
%             % process
%             c = regexpcellout(b,searchterm,'match');
%             i = regexpcellout(b,searchterm);
% %                     b(i) = [];
%             c(cellfun(@isempty,c)) = [];
%             if numel(c)>1
%                 c = strjoinrows(c',', ');
%             end
% 
%             if ~isempty(c)
%                 T.(name)(ai) = c;
%             end
%         end
% return
%     end
% end


%% go through each strain
% varnames = {'genotype','gene','promoter','allele','chromosome','array','others'};
% A = cell(size(genotypes,1),numel(varnames));
% T = cell2table(A,'VariableNames',varnames);
% startn = 1;
% T.genotype = genotypes;
% 
% 
% for vi = 2:numel(varnames)
%     name = varnames{vi};
%     
%     switch name
%         case 'gene'
%             a = regexpcellout(genotypes,'[;]','split');
%             % extract genes
%             for ai = 1:size(a,2)
%                 b = a(:,ai);
%                 % delete allele (
%                 c = regexpcellout(b,'[(]','split');
%                 b = c(:,1);
%                 % delete before ::
%                 b = regexprep(b,'.+(?=[:])','');
%                 % delete : and ]
%                 b = regexprep(b,'(:)|(])','');
%                 % pull out only gene info
%                 searchterm = '([a-zA-Z]{1,}-\d+)|(^wildtype$)|(^[A-Z]\d+[A-Z]+\d+[.]\d+)';
%                 c = regexpcellout(b,searchterm,'match');
%                 % put back in a
%                 if size(c,2) >1
%                     i = find(~cellfun(@isempty,c(:,2)));
%                     for ii = 1:numel(i)
%                         d = c(i(ii),:);
%                         d(cellfun(@isempty,d)) = [];
%                         c(i) = strjoinrows(d ,', ');
%                     end
%                 end
%                 a(:,ai) = c(:,1);
%             end
%             % remove duplicate and join
%            for ai = 1:size(a,1)
%                b = a(ai,:);
%                b(cellfun(@isempty,b)) = [];
%                b = unique(b);
%                a(ai,1) = strjoinrows(b,', ');
%            end
%            o = a(:,1);
% 
%         case 'chromosome'
%             a = regexpcellout(genotypes,'[(]|[)]|[;]|[[]|[]]','split');
%             a = genotypes;
%             searchterm = '(?<=[)])((I+)|IV|X|V)';
%             a = regexpcellout(a,searchterm,'match');
%             % remove duplicate and join
%            for ai = 1:size(a,1)
%                b = a(ai,:);
%                b(cellfun(@isempty,b)) = [];
%                a(ai,1) = strjoinrows(b,', ');
%            end
%            o = a(:,1);
% 
%         case 'promoter'
%             a = regexpcellout(genotypes,'[(]|[)]|[;]|[[]|[]]','split');
%             searchterm = '(?<=[;])(P|p)[a-zA-Z]+[-]\d+(??=[:])';
%             a = genotypes;
%             searchterm = '(?<=[)])((I+)|IV|X|V)';
%             a = regexpcellout(a,searchterm,'match');
%             % remove duplicate and join
%            for ai = 1:size(a,1)
%                b = a(ai,:);
%                b(cellfun(@isempty,b)) = [];
%                a(ai,1) = strjoinrows(b,', ');
%            end
%            o = a(:,1);
%            
%         case 'array'
%             searchterm = '(^(EX)|(Ex)([a-z]{1,}|[0-9]{1,}))|([a-z]{1,}Is[0-9]{1,})';
%         case 'allele'
%             searchterm = '^[a-z]{1,}\d{1,}|(hawaiian)|(bristol)';
%     end
% %     switch name 
% %         case 'others'
% %             c = strjoinrows(b',', ');
% %         otherwise
% %             % process
% %             c = regexpcellout(b,searchterm,'match');
% %             i = regexpcellout(b,searchterm);
% % %             b(i) = [];
% % 
% %             c(cellfun(@isempty,c)) = [];
% %             if numel(c)>1
% %                 c = strjoinrows(c',', ');
% %             end
% % 
% %     end
% %     if ~isempty(c)
%         T.(name) = o;
% %     end
% 

% end




%% reorganize table
% T1 = innerjoin(S,T);
% writetable(T1,fullfile(pM,'strain information.csv'));