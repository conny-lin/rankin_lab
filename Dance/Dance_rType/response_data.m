function [baseline, response, response_lower, response_rev, response_higher] = response_data(GraphData,fn, xtime, rtime, baselinetime)

%%
    baseline = GraphData.(fn)(xtime == baselinetime);
    response = GraphData.(fn)(ismember(xtime,rtime));
    response_lower = response < baseline;
    response_rev = response < 0;
    response_higher = response > baseline;

    