function MTemp = MWTDBquery_grpInSameExp(MWTDB,rgroup,cgroup,varargin)


    searchkey = 'groupname';
    vararginProcessor;
%     rgroup = {strain,[strain,'_400mM']};
%     cgroup = {'N2','N2_400mM'};
%     
    i = ismember(MWTDB.(searchkey),rgroup);
    i = ismember(MWTDB.expname, MWTDB.expname(i));
    MTemp = MWTDB(i,:);
    MTemp(~ismember(MTemp.(searchkey),[rgroup;cgroup]),:) = [];
    
    
end