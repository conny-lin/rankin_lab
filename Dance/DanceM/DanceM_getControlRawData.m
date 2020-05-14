function [ControlData,pMWTC] = DanceM_getControlRawData(RawData,pMWT,ctrl_groupname)


MWTDB = parseMWTinfo(pMWT);
groupname = MWTDB.groupname;

%% calculate control mean by experiment
i_ctrl = ismember(groupname,ctrl_groupname);
pMWTC = pMWT(i_ctrl);
% get msr list
msrlist = fieldnames(RawData);

% get control data
ControlData = struct;
for msri = 1:numel(msrlist)
    msrname = msrlist{msri};
    ControlData.(msrname) = RawData.(msrname)(:,i_ctrl);
end