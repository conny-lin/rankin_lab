function MWTDB = parseToleranceName(MWTDB)


% cal
T = parseRx(MWTDB.rx);

mM = repmat({'mM'},size(T,1),1); % mM text
predose = strjoinrows([num2cellstr(T.dose_rx) mM],'');
postdose = strjoinrows([num2cellstr(T.dose_test) mM],'');
cond = strjoinrows([predose postdose]);
MWTDB.groupname_short = strjoinrows([MWTDB.strain cond]);
MWTDB.condition_short = cond;
MWTDB.predose = predose;
MWTDB.postdose = postdose;