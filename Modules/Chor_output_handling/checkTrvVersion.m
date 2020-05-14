function trv_oldversion = checkTrvVersion(pTrv)

trv_oldversion = true(size(pTrv));
for mwti =1:numel(pTrv)
    pf = pTrv{mwti};
    % see version of trv
    fileID = fopen(pf,'r');
    a = textscan(fileID,'%s', 1,'Delimiter', '', 'WhiteSpace', '');
    a = char(a{1}(1));
    fc = a(1);
    fclose(fileID);
    if regexp(fc,'\d{1,}') == 1; 
        trv_oldversion(mwti) = false;
    end
end