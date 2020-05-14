function [S,A] = cal_velocity(S)



%% calculate velocity
Bias = S.bias.array;
Speed = S.speed.array;
A = (Speed.*Bias); % calculate output
% fill in velocity for bias == 0
% A(Bias==0) = Speed(Bias==0); % this one can cause issue

%%
S.velocity.array = A;
S.velocity.array_transposed = A.';
S.velocity.timetable = array2timetable(A,'RowTimes',S.bias.timetable.Time,'VariableNames',S.bias.timetable.Properties.VariableNames);
A = S.velocity.array;

