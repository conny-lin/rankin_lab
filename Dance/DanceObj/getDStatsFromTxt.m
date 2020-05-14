function T = getDStatsFromTxt(c)

    % get data from text file +++++++++++++++
    T = table;
    % get group name
    a = c(regexpcellout(c,'gname = '));
    b = regexpcellout(a,'(=)|(,)','split');
    b = b(:,2:end);
    b = regexprep(b,' ','_');
    b = regexprep(b,'^[_]','');
    T.gname = b';
    % get stats
    stnames = {'n','mean','se'};
    for x = 1:numel(stnames)
        % get n
        stname = stnames{x};
        a = c(regexpcellout(c,sprintf('^%s = ',stname)));
        b = regexpcellout(a,'(=)|(,)','split');
        b = b(:,2:end);
        b = regexprep(b,' ','');
        T.(stname) = cellfun(@str2num,b)';
    end