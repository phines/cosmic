clear all; close all; 

% load nminus2set;
outputs = [];
psout   = [];
m       = 2383;

% for i = 1:size(BOpairs,1)
%     disp(i);
%     a = BOpairs(i,1);
%     b = BOpairs(i,2);
%     [outputs{i},psout{i}] = sim_case2383_n_2(a,b);
% end
% save results_case2383_nminus2 outputs psout;

for i = 1:m
    a = i;
    fprintf('running contingency #%d: branches %d ...\n',i,a);
    [outputs{i},psout{i}] = sim_case2383_n_2(a);
end
save results_case2383_nminus1 outputs psout;