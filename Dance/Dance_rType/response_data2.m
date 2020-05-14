function [rqualifier, rdata] = response_data2(GraphData,gname)

%%
alpha = 0.05;
pvlimit = 0.001;
responseT1 = 0.1;
responseT2 = 1;
baselinetime = -0.1;
xtime = GraphData.N2_x;
rtime = GraphData.N2_x(GraphData.N2_x >= responseT1 & ...
        GraphData.N2_x <= responseT2);
    
fn = [gname,'_y'];

%%
rdata = struct;
rdata.baseline = GraphData.(fn)(xtime == baselinetime);
rdata.response = GraphData.(fn)(ismember(xtime,rtime));
rqualifier.lower = rdata.response < rdata.baseline;
rqualifier.rev = rdata.response < 0;
rqualifier.higher = rdata.response > rdata.baseline;

    