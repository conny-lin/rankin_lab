function [preplate,ISI,tapN,postrecord] = parse_runcond(runcond)

%% run condition examples  
% runcond = '100s30x10s10s';
% runcond = {runcond;runcond}


%% split
if iscell(runcond) == 0
    runcond = {runcond};
end

a = regexpcellout(runcond,'s|x','split');
preplate = str2double(a{:,1});
ISI = str2double(a{:,3});
tapN = str2double(a{:,2});
postrecord = str2double(a{:,4});


