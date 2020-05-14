function MWTDB = rxname_add2MWTDB(MWTDB)

% change to group name alternate for igor
cd('/Users/connylin/Dropbox/Code/Matlab/Library RL/Modules/MWTDatabase')
rxname = readtable('rxname_alternate.csv');
for ri =1:size(rxname,1)
    i= ismember(MWTDB.rx,rxname.rx(ri));
    MWTDB.rx(i) = rxname.rx_alternate(ri);
end
% create new names
A = [MWTDB.strain MWTDB.rx];
MWTDB.groupname = strjoinrows(A,'_');
MWTDB = [MWTDB parseRx(MWTDB.rx)]; % add treatment conditions
