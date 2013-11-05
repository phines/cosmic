clear all; close all; rng(28);

n_reps  = 100;
m       = 46;
outputs = [];
psout   = [];

for i = 1:n_reps
    a = randi(m,1); b = a;
    while b == a
        b = randi(m,1);
    end
    fprintf('running contingency #%d: branches %d and %d ...\n',i,a,b);
    [outputs{i},psout{i}] = sim_case39_n_2(a,b);
end

save results_case39_nminus2 outputs psout;
