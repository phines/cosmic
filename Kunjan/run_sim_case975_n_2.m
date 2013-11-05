clear all; close all; 

n_reps  = 100;
m       = 784;
outputs = [];
psout   = [];

for i = 1:n_reps
    a = randi(m,1); b = a;
    while b == a
        b = randi(m,1);
    end
    fprintf('running contingency #%d: branches %d and %d ...\n',i,a,b);
    [outputs{i},psout{i}] = sim_case975_n_2(a,b);
end

save results_case975_nminus2 outputs psout;
