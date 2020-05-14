function str = sample_size_text(strain)

pSampleSize = '/Users/connylin/Dropbox/RL Pub PhD Dissertation/Chapters/4-EARS/3-Results/2-Genes EARS/Summary v3/ephys_accpeak_samplesize';


%% sample size
p = fullfile(pSampleSize, [strain, ' ephys t28-30 samplesize.csv']);
T = readtable(p);
n_sample = T.wormN;
nstring = strjoin(num2cellstr(n_sample)',', ');

% sample name
sample_name = T.gname;
% add 0mM
i = ~regexpcellout(sample_name,'400mM');
a = sample_name(i);
b = repmat({'0mM'},size(a));
c = strjoinrows([a,b],' ');
sample_name(i) = c;
% replace _
sample_name = regexprep(sample_name,'_',' ');
sample_name_str = strjoin(sample_name',', ');

% get n strain
str = sprintf('N(%s)=%s',sample_name_str ,nstring);
