clear all; close all; rng(28);

n_reps  = 300;
m       = 46;
% n_reps  = combnk(1:m,2);
outputs = [];
psout   = [];

% for i = 1:n_reps
% % for i = 1:size(n_reps,1)
% %     comb = n_reps(i,:);
%     a = randi(m,1); b = a;
%     while b == a
%         b = randi(m,1);
%     end
% %     a = comb(1); b = comb(2);
%     fprintf('running contingency #%d: branches %d and %d ...\n',i,a,b);
%     [outputs{i},psout{i}] = sim_case39_n_2(a,b);
% end
% 
% save results_case39_nminus2_polar outputs psout;

for i = 1:m
    a = i;
    fprintf('running contingency #%d: branches %d ...\n',i,a);
    [outputs{i},psout{i}] = sim_case39_n_2(a);
end
save results_case39_nminus1_polar outputs psout;
% save results_case39_nminus1_rec outputs psout;
