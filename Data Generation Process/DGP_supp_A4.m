global  n x m Rep

Rep = 500; % # of Monte Carlo replication
n_choice = [120; 240];

p = size(n_choice,1);

seed = 301;
beta0 = 0;
sigma = 1;

for i = 1:p
    n = n_choice(i);
    for r = 1:Rep
        seed_v = seed + r;
        rng(seed_v);
        x = random('Normal', beta0, sigma, n, 1);
        eval(['data_n_',num2str(n), '(',num2str(r),',','1' , ')=','{x}']);
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global  n x m Rep

Rep = 500; % # of Monte Carlo replication
n_choice = [120; 240];

p = size(n_choice,1);

seed = 301;
theta0 = 0.5;

for i = 1:p
    n = n_choice(i);
    for r = 1:Rep
        seed_v = seed + r;
        rng(seed_v);
        E = random('Normal', 0, 0.9, n, 1);
        Z = random('Normal', 0, 1, n, 1);
        y = Z*theta0 + E;
        eval(['data_n_',num2str(n), '(',num2str(r),',','1' , ')=','{[y, Z]}']); 
    end
end





