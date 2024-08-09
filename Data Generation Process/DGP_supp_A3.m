global  n y x z m Rep

Rep = 500; % # of Monte Carlo replication
n_choice = [20; 40; 60; 80; 100; 120];
m_choice = [80; 160];
p = size(n_choice,1);
q = size(m_choice,1);

seed = 301;
useful = 4;

rho = 0.6; % the endogeneity
beta0 = [0.5; 0.5];

for i = 1:p
    n = n_choice(i);
    for j = 1:q
        m = m_choice(j);
        eval(['data_n_', num2str(n),'_m_', num2str(m) '=', 'cell(500,1)' ]);
        for r = 1:Rep
            %Fix the seed for each iteration
            disp(r);
            seed_v = seed + r;
            rng(seed_v);
            [y, x, z] = dgpLinearIV(beta0, rho, useful); % generate the data
            eval(['data_n_',num2str(n),'_m_', num2str(m), '(',num2str(r),',','1' , ')=','{[y, x, z]}']);
        end
    end
end


