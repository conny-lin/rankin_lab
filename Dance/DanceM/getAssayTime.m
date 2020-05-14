function TimeSet = getAssayTime(Data,assaytype,varargin)


%% get time assay information
TimeSet = struct;
TimeSet.assaytype = assaytype;


switch assaytype
    case 'rType'
        % input data is MWTDB
        T = Data(:,{'ISI','preplate','postrec','tapN'});
        TU = unique(T);
        if size(TU,1)~=1
            disp(TU);
            error('can not handle more than one exp set up');
        end
        
        tf = ((TU.tapN-1)*TU.ISI)+TU.preplate;
        taptimes = TU.preplate : TU.ISI : tf;
        beforetap = 0.3;
        aftertap = 1;
        atstart = taptimes-beforetap;
        atend = taptimes+aftertap;
        
        TimeSet.taptimes = taptimes;
        TimeSet.atstart = atstart;
        TimeSet.atend = atend;
        TimeSet.beforetap = beforetap;
        TimeSet.aftertap = aftertap;

end


