clear all
close all
clc

%% Price Sensitivity
tic;

% Initilization
interval_num = 2000;
sample_size = 1e6;

% Price Sensitivity Distribution
nd = makedist('Normal', 'mu', 1, 'sigma', 0.5);
price_sensitivity_dist = truncate(nd, 0.5, 1.5);
x = linspace(0.5,1.5,interval_num);
res = pdf(price_sensitivity_dist, x);
price_sensitivity_rand = randpdf(res, x, [sample_size,1]);

% Plot
fig = figure;
histogram(price_sensitivity_rand,50);
title('Price Sensitivity');

% Save Figure
saveas(fig,'./Figures/price_sensitivity_sim1.png');

% end
toc;

%% Consumption Preference
tic;

% Initilization
interval_num = 2000;
sample_size = 1e6;

% Consumption Preference Distribution
nd = makedist('Normal', 'mu', 30, 'sigma', 4);
cons_pref_dist = truncate(nd, 0, 60);
x = linspace(15,45,interval_num);
res = pdf(cons_pref_dist, x);
cons_pref_rand = randpdf(res, x, [sample_size,1]);

% Plot
fig = figure;
histogram(cons_pref_rand,50);
title('Consumption Preference');

% Save Figure
saveas(fig,'./Figures/consumption_preference_sim1.png');

% end
toc;

%% Clear Variables - Stop 1
clear fig interval_num nd res x sample_size

%% Simulation - Consumption per Capita
tic;

% Initilization
Hours = 6;
run_time = 50;
N = 100;

% Hybrid Contents
hybrid = [0.8, 0.1, 0.1, 0; 0.1, 0.8, 0.1, 0; 0.1, 0.1, 0.8, 0; 0, 0.1, 0.1, 0.8;...
    0, 0.8, 0.1, 0.1; 0, 0.1, 0.8, 0.1; 0.20, 0.30, 0.15, 0.35; 0.15, 0.35, 0.20, 0.30;...
    0.25, 0.35, 0.15, 0.35; 0.20, 0.30, 0.20, 0.30];
hybrid = hybrid * N;

% Total Consumption Matricies
total_consumption_S = zeros(3,run_time,Hours);
total_consumption_U = zeros(3,run_time,Hours);
total_consumption_W = zeros(3,run_time,Hours);
total_consumption_C = zeros(3,run_time,Hours);
total_consumption_H1 = zeros(3,run_time,Hours);
total_consumption_H2 = zeros(3,run_time,Hours);
total_consumption_H3 = zeros(3,run_time,Hours);
total_consumption_H4 = zeros(3,run_time,Hours);
total_consumption_H5 = zeros(3,run_time,Hours);
total_consumption_H6 = zeros(3,run_time,Hours);
total_consumption_H7 = zeros(3,run_time,Hours);
total_consumption_H8 = zeros(3,run_time,Hours);
total_consumption_H9 = zeros(3,run_time,Hours);
total_consumption_H10 = zeros(3,run_time,Hours);

% Estimated and Realized Value Matricies
estimated_total_consumption_container = zeros(run_time,Hours);
estimated_total_revenue_container = zeros(run_time,Hours);
estimated_total_cost_container = zeros(run_time,Hours);
estimated_profit_container = zeros(run_time,1);
estimated_welfare_container = zeros(run_time,Hours);

realized_profit_rate = zeros(run_time,42,Hours);
price_container = zeros(run_time,42,Hours,N);

% Cost Containers
cost_container_S = zeros(run_time,3,Hours,N);
cost_container_U = zeros(run_time,3,Hours,N);
cost_container_W = zeros(run_time,3,Hours,N);
cost_container_C = zeros(run_time,3,Hours,N);
cost_container_H1 = zeros(run_time,3,Hours,N);
cost_container_H2 = zeros(run_time,3,Hours,N);
cost_container_H3 = zeros(run_time,3,Hours,N);
cost_container_H4 = zeros(run_time,3,Hours,N);
cost_container_H5 = zeros(run_time,3,Hours,N);
cost_container_H6 = zeros(run_time,3,Hours,N);
cost_container_H7 = zeros(run_time,3,Hours,N);
cost_container_H8 = zeros(run_time,3,Hours,N);
cost_container_H9 = zeros(run_time,3,Hours,N);
cost_container_H10 = zeros(run_time,3,Hours,N);

% Revenue Containers
revenue_container_S = zeros(run_time,3,Hours);
revenue_container_U = zeros(run_time,3,Hours);
revenue_container_W = zeros(run_time,3,Hours);
revenue_container_C = zeros(run_time,3,Hours);
revenue_container_H1 = zeros(run_time,3,Hours);
revenue_container_H2 = zeros(run_time,3,Hours);
revenue_container_H3 = zeros(run_time,3,Hours);
revenue_container_H4 = zeros(run_time,3,Hours);
revenue_container_H5 = zeros(run_time,3,Hours);
revenue_container_H6 = zeros(run_time,3,Hours);
revenue_container_H7 = zeros(run_time,3,Hours);
revenue_container_H8 = zeros(run_time,3,Hours);
revenue_container_H9 = zeros(run_time,3,Hours);
revenue_container_H10 = zeros(run_time,3,Hours);

% Welfare Containers
welfare_S = zeros(run_time,3,Hours);
welfare_U = zeros(run_time,3,Hours);
welfare_W = zeros(run_time,3,Hours);
welfare_C = zeros(run_time,3,Hours);
welfare_H1 = zeros(run_time,3,Hours);
welfare_H2 = zeros(run_time,3,Hours);
welfare_H3 = zeros(run_time,3,Hours);
welfare_H4 = zeros(run_time,3,Hours);
welfare_H5 = zeros(run_time,3,Hours);
welfare_H6 = zeros(run_time,3,Hours);
welfare_H7 = zeros(run_time,3,Hours);
welfare_H8 = zeros(run_time,3,Hours);
welfare_H9 = zeros(run_time,3,Hours);
welfare_H10 = zeros(run_time,3,Hours);

% Consumption
consumption_S = zeros(run_time,3,N,Hours);
consumption_U = zeros(run_time,3,N,Hours);
consumption_W = zeros(run_time,3,N,Hours);
consumption_C = zeros(run_time,3,N,Hours);
consumption_H1 = zeros(run_time,3,N,Hours);
consumption_H2 = zeros(run_time,3,N,Hours);
consumption_H3 = zeros(run_time,3,N,Hours);
consumption_H4 = zeros(run_time,3,N,Hours);
consumption_H5 = zeros(run_time,3,N,Hours);
consumption_H6 = zeros(run_time,3,N,Hours);
consumption_H7 = zeros(run_time,3,N,Hours);
consumption_H8 = zeros(run_time,3,N,Hours);
consumption_H9 = zeros(run_time,3,N,Hours);
consumption_H10 = zeros(run_time,3,N,Hours);

% Utility Container
utility_S = zeros(run_time,3,N,Hours);
utility_U = zeros(run_time,3,N,Hours);
utility_W = zeros(run_time,3,N,Hours);
utility_C = zeros(run_time,3,N,Hours);
utility_H1 = zeros(run_time,3,N,Hours);
utility_H2 = zeros(run_time,3,N,Hours);
utility_H3 = zeros(run_time,3,N,Hours);
utility_H4 = zeros(run_time,3,N,Hours);
utility_H5 = zeros(run_time,3,N,Hours);
utility_H6 = zeros(run_time,3,N,Hours);
utility_H7 = zeros(run_time,3,N,Hours);
utility_H8 = zeros(run_time,3,N,Hours);
utility_H9 = zeros(run_time,3,N,Hours);
utility_H10 = zeros(run_time,3,N,Hours);

% Other Matricies
estimated_revenue = zeros(run_time,Hours);
estimated_cost = zeros(run_time,Hours);
estimated_consumption = zeros(run_time,N,Hours);
estimated_total_consumption = zeros(run_time,Hours);
expected_total_consumption_squared = zeros(run_time,Hours);

% Other Parameters and Configurations
realized_omega = zeros(run_time,1);
seed= rng;
rng('default')
sigma_pop = 1;
sigma_ij = 0;

for run = 1:run_time
    current_sigma_ij = sigma_ij;
    time_zone = ones(1,Hours);
    
    graph = generate_graph(N,run);
    incidence_matrix = graph;
    nr_nodes = N;
    nr_neighbors = incidence_matrix*ones(nr_nodes,1);
    
    beta = zeros(Hours,1);
    
    % Setup
    g_bar = mean(cons_pref_rand);
    g = zeros(1,N);
    e = zeros(1,N);
    sigma_ii = 4;
    Sigma = sigma_ii * eye(N,N) + current_sigma_ij * (ones(N,N)-eye(N,N));
    
    % Inference Matrix
    IM = zeros(N,N);
    for i = 1:N
        for j = 1:N
            if i ~= j
                IM(i,j) = Sigma(i,j)/Sigma(i,i);
            end
        end
    end
    
    omega_bar = 0;
    omega_variance = 2;
    profit_rate = 1.2; %Fixed profit rate
    kappa = 1;
    beta(:) = profit_rate * kappa * ones(Hours,1);
    beta_container(:,run) = beta(:);
    sigma_ij = current_sigma_ij/sigma_ii;
    
    % Estimation of consumption, cost etc.
    for h = 1:Hours
       k = 0;
       if time_zone(1,h) == 1
           k = 1;
       end

       a = ones(N,1)/((N+1)*beta(h) + 2*mean(cons_pref_rand));
       rho = 1/(2*(beta(h) + mean(cons_pref_rand)));
       b = ones(N,1)*rho/((N-1)*beta(h)*rho*sigma_ij+1);

       estimated_consumption(run,:,h) = a*(g_bar(k)-omega_bar*beta(h));
       estimated_total_consumption(run,h) = sum(estimated_consumption(run,:,h));
       expected_total_consumption_squared(run,h) = estimated_total_consumption(run,h)...
           *estimated_total_consumption(run,h) + sigma_ii*sigma_ii* b'*b ;
       estimated_revenue(run,h) = beta(h) * expected_total_consumption_squared(run,h) ...
           + omega_bar *estimated_total_consumption(run,h);
       estimated_cost(run,h) = kappa *expected_total_consumption_squared(run,h);
       expected_profit_rate(run,1) = estimated_revenue(run,h)/estimated_cost(run,h);
    end
    
    estimated_total_consumption_container(run,:) =  estimated_total_consumption(run,:);
    estimated_total_revenue_container(run,:) = estimated_revenue(run,:);
    estimated_total_cost_container(run,:) = estimated_cost(run,:);
    estimated_profit_container(run) = expected_profit_rate(run);
    total_estimated_cost = sum(estimated_cost);
    
    % Simulation
    rng(run)
    rand_indices = ceil(length(cons_pref_rand)*rand(1,N));
    g(1,:) = cons_pref_rand(rand_indices);
    e(1,:) = price_sensitivity_rand(rand_indices);
    
    omega = omega_bar + omega_variance * randn(1);
    realized_omega(run) = omega;
    cost = zeros(Hours,1);
    revenue = zeros(Hours,1);
    price = zeros(Hours,N);
    
    for h = 1:Hours
       k=0;
       if time_zone(1,h) == 1
           k = 1;
       end
       
       % Selfish - Private
       tau = beta(h)/N; lambda = beta(h)/N; mu = beta(h)/N; %main parameters

       rho = 1./(2.*(tau+e));
       a_SP = rho./(1 + lambda.*rho.*sigma_ij.*(N-1));
       b_SP = rho./(1 + lambda.*rho.*(N-1));

       consumption_S(run,1,:,h) = b_SP'.*(g_bar(k)-omega_bar*mu) + a_SP'.*(g(k,:) - g_bar(k))';
       total_consumption_S(1,run,h) = sum(consumption_S(run,1,:,h));
       price(h,:) = beta(h) * (total_consumption_S(1,run,h) + omega)/N;
       price_container(run,1,h,:) = price(h,:);

       cost(h) = kappa * total_consumption_S(1,run,h) ...
           * total_consumption_S(1,run,h)/N;
       revenue(h) = price(h,1)* total_consumption_S(1,run,h);
       realized_profit_rate(run,1,h) = revenue(h)/ cost(h);
       cost_container_S(run,1,h,1) = cost(h);

       revenue_container_S(run,1,h) = revenue(h);
       for i = 1:N
           u_i = -consumption_S(run,1,i,h)*(beta(h)*(total_consumption_S(1,run,h)...
               +omega)/N)+g(k,i)*consumption_S(run,1,i,h)-e(i)*consumption_S(run,1,i,h)*consumption_S(run,1,i,h);
           utility_S(run,1,i,h) = u_i;
           welfare_S(run,1,h) = welfare_S(run,1,h) + u_i;
       end

       % Selfish - Broadcast
       tau = beta(h)/N; lambda = beta(h)/N; mu = beta(h)/N; %main parameters

       rho = 1./(2.*(tau+e));
       a_SB = rho./(1 + lambda.*rho);
       b_SB = rho./(1 + lambda.*rho.*(N-1));

       consumption_S(run,2,:,h) = a_SB.*(g(k,:)-sum(g(k,:))/N) + b_SB.*(sum(g(k,:))/N-omega_bar*mu);
       total_consumption_S(2,run,h) = sum(consumption_S(run,2,:,h));
       price(h,:) = beta(h) * (total_consumption_S(2,run,h) + omega)/N;
       price_container(run,2,h,:) = price(h,:);

       cost(h) = kappa * total_consumption_S(2,run,h) ...
           * total_consumption_S(2,run,h)/N;
       revenue(h) = price(h,1)* total_consumption_S(2,run,h);
       realized_profit_rate(run,2,h) = revenue(h)/ cost(h);
       cost_container_S(run,2,h,1) = cost(h);

       revenue_container_S(run,2,h) = revenue(h);
       for i = 1:N
           u_i = -consumption_S(run,2,i,h)*(beta(h)*(total_consumption_S(2,run,h)...
               +omega)/N)+g(k,i)*consumption_S(run,2,i,h)-e(i)*consumption_S(run,2,i,h)*consumption_S(run,2,i,h);
           utility_S(run,2,i,h) = u_i;
           welfare_S(run,2,h) = welfare_S(run,2,h) + u_i;
       end
       
       % Altruistic - Private
       tau = beta(h)/N; lambda = 2*beta(h)/N; mu = beta(h)/N; %main parameters

       rho = 1./(2.*(tau+e));       
       a_UP = rho./(1 + lambda.*rho.*sigma_ij.*(N-1));
       b_UP = rho./(1 + lambda.*rho.*(N-1));

       consumption_U(run,1,:,h) = b_UP'.*(g_bar(k)-omega_bar*mu) + a_UP'.*(g(k,:) - g_bar(k))';
       total_consumption_U(1,run,h) = sum(consumption_U(run,1,:,h));
       price(h,:) = beta(h) * (total_consumption_U(1,run,h) + omega)/N;
       price_container(run,4,h,:) = price(h,:);

       cost(h) = kappa * total_consumption_U(1,run,h) ...
           * total_consumption_U(1,run,h)/N;
       revenue(h) = price(h,1)* total_consumption_U(1,run,h);
       realized_profit_rate(run,4,h) = revenue(h)/ cost(h);
       cost_container_U(run,1,h,1) = cost(h);

       revenue_container_U(run,1,h) = revenue(h);
       for i = 1:N
           u_i = -consumption_U(run,1,i,h)*(beta(h)*(total_consumption_U(1,run,h)...
               +omega)/N)+g(k,i)*consumption_U(run,1,i,h)-e(i)*consumption_U(run,1,i,h)*consumption_U(run,1,i,h);
           utility_U(run,1,i,h) = u_i;
           welfare_U(run,1,h) = welfare_U(run,1,h) + u_i;
       end

       % Altruistic - Broadcast
       tau = beta(h)/N; lambda = 2*beta(h)/N; mu = beta(h)/N; %main parameters

       rho = 1./(2.*(tau+e));
       a_UB = rho./(1-lambda.*rho);
       b_UB = rho./(1+lambda.*rho.*(N-1));

       consumption_U(run,2,:,h) = a_UB.*(g(k,:)-sum(g(k,:))/N) + b_UB.*(sum(g(k,:))/N-omega_bar*mu);
       total_consumption_U(2,run,h) = sum(consumption_U(run,2,:,h));
       price(h,:) = beta(h) * (total_consumption_U(2,run,h) + omega)/N;
       price_container(run,5,h,:) = price(h,:);

       cost(h) = kappa * total_consumption_U(2,run,h) ...
           * total_consumption_U(2,run,h)/N;
       revenue(h) = price(h,1)* total_consumption_U(2,run,h);
       realized_profit_rate(run,5,h) = revenue(h)/ cost(h);
       cost_container_U(run,2,h,1) = cost(h);

       revenue_container_U(run,2,h) = revenue(h);
       for i = 1:N
           u_i = -consumption_U(run,2,i,h)*(beta(h)*(total_consumption_U(2,run,h)...
               +omega)/N)+g(k,i)*consumption_U(run,2,i,h)-e(i)*consumption_U(run,2,i,h)*consumption_U(run,2,i,h);
           utility_U(run,2,i,h) = u_i;
           welfare_U(run,2,h) = welfare_U(run,2,h) + u_i;
       end
       
       % Welfare Maximizer - Private
       tau = kappa/N; lambda = 2*kappa/N; mu = 0; %main parameters

       rho = 1./(2.*(tau+e));       
       a_WP = rho./(1 + lambda.*rho.*sigma_ij.*(N-1));
       b_WP = rho./(1 + lambda.*rho.*(N-1));

       consumption_W(run,1,:,h) = b_WP'.*(g_bar(k)-omega_bar*mu) + a_WP'.*(g(k,:) - g_bar(k))';
       total_consumption_W(1,run,h) = sum(consumption_W(run,1,:,h));
       price(h,:) = beta(h) * (total_consumption_W(1,run,h) + omega)/N;
       price_container(run,7,h,:) = price(h,:);

       cost(h) = kappa * total_consumption_W(1,run,h) ...
           * total_consumption_W(1,run,h)/N;
       revenue(h) = price(h,1)* total_consumption_W(1,run,h);
       realized_profit_rate(run,7,h) = revenue(h)/ cost(h);
       cost_container_W(run,1,h,1) = cost(h);

       revenue_container_W(run,1,h) = revenue(h);
       for i = 1:N
           u_i = -consumption_W(run,1,i,h)*(beta(h)*(total_consumption_W(1,run,h)...
               +omega)/N)+g(k,i)*consumption_W(run,1,i,h)-e(i)*consumption_W(run,1,i,h)*consumption_W(run,1,i,h);
           utility_W(run,1,i,h) = u_i;
           welfare_W(run,1,h) = welfare_W(run,1,h) + u_i;
       end

       % Welfare Maximizer - Broadcast
       tau = kappa/N; lambda = 2*kappa/N; mu = 0; %main parameters

       rho = 1./(2.*(tau+e));
       a_WB = rho./(1 - lambda.*rho);
       b_WB = rho./(1 + lambda.*rho.*(N-1));

       consumption_W(run,2,:,h) = a_WB.*(g(k,:)-sum(g(k,:))/N) + b_WB.*(sum(g(k,:))/N-omega_bar*mu);
       total_consumption_W(2,run,h) = sum(consumption_W(run,2,:,h));
       price(h,:) = beta(h) * (total_consumption_W(2,run,h) + omega)/N;
       price_container(run,8,h,:) = price(h,:);

       cost(h) = kappa * total_consumption_W(2,run,h) ...
           * total_consumption_W(2,run,h)/N;
       revenue(h) = price(h,1)* total_consumption_W(2,run,h);
       realized_profit_rate(run,8,h) = revenue(h)/ cost(h);
       cost_container_W(run,2,h,1) = cost(h);

       revenue_container_W(run,2,h) = revenue(h);
       for i = 1:N
           u_i = -consumption_W(run,2,i,h)*(beta(h)*(total_consumption_W(2,run,h)...
               +omega)/N)+g(k,i)*consumption_W(run,2,i,h)-e(i)*consumption_W(run,2,i,h)*consumption_W(run,2,i,h);
           utility_W(run,2,i,h) = u_i;
           welfare_W(run,2,h) = welfare_W(run,2,h) + u_i;
       end
       
       % Community Friendly - Private
       tau = beta(h)/N; mu = beta(h)/N; lambda = (N + nr_neighbors')*beta(h)/(N^2); %main parameters
       rho = 1./(2.*(tau+e));       
       a_CP = rho./(1 + lambda.*rho.*sigma_ij.*(N-1));
       b_CP = rho./(1 + lambda.*rho.*(N-1));

       consumption_C(run,1,:,h) = b_CP'.*(g_bar(k)-omega_bar*mu) + a_CP'.*(g(k,:) - g_bar(k))';
       total_consumption_C(1,run,h) = sum(consumption_C(run,1,:,h));
       price(h,:) = beta(h) * (total_consumption_C(1,run,h) + omega)/N;
       price_container(run,10,h,:) = price(h,:);

       cost(h) = kappa * total_consumption_C(1,run,h) ...
           * total_consumption_C(1,run,h)/N;
       revenue(h) = price(h,1)* total_consumption_C(1,run,h);
       realized_profit_rate(run,10,h) = revenue(h)/ cost(h);
       cost_container_C(run,1,h,1) = cost(h);

       revenue_container_C(run,1,h) = revenue(h);
       for i = 1:N
           u_i = -consumption_C(run,1,i,h)*(beta(h)*(total_consumption_C(1,run,h)...
               +omega)/N)+g(k,i)*consumption_C(run,1,i,h)-e(i)*consumption_C(run,1,i,h)*consumption_C(run,1,i,h);
           utility_C(run,1,i,h) = u_i;
           welfare_C(run,1,h) = welfare_C(run,1,h) + u_i;
       end

       % Community Friendly - Broadcast
       tau = beta(h)/N; mu = beta(h)/N; lambda = (N + nr_neighbors')*beta(h)/(N^2); %main parameters
       rho = 1./(2.*(tau+e));
       a_CB = rho./(1 - lambda.*rho);
       b_CB = rho./(1 + lambda.*rho.*(N-1));

       consumption_C(run,2,:,h) = a_CB.*(g(k,:)-sum(g(k,:))/N) + b_CB.*(sum(g(k,:))/N-omega_bar*mu);
       total_consumption_C(2,run,h) = sum(consumption_C(run,2,:,h));
       price(h,:) = beta(h) * (total_consumption_C(2,run,h) + omega)/N;
       price_container(run,11,h,:) = price(h,:);

       cost(h) = kappa * total_consumption_C(2,run,h) ...
           * total_consumption_C(2,run,h)/N;
       revenue(h) = price(h,1)* total_consumption_C(2,run,h);
       realized_profit_rate(run,11,h) = revenue(h)/ cost(h);
       cost_container_C(run,2,h,1) = cost(h);

       revenue_container_C(run,2,h) = revenue(h);
       for i = 1:N
           u_i = -consumption_C(run,2,i,h)*(beta(h)*(total_consumption_C(2,run,h)...
               +omega)/N)+g(k,i)*consumption_C(run,2,i,h)-e(i)*consumption_C(run,2,i,h)*consumption_C(run,2,i,h);
           utility_C(run,2,i,h) = u_i;
           welfare_C(run,2,h) = welfare_C(run,2,h) + u_i;
       end    
    end
    
    for hyb = 1:10
        customers = zeros(1,N);
        sp = 1;
        for nums = 1:4
            if sp > 100
                break
            end
            customers(sp:sp+hybrid(hyb,nums)-1) = nums;
            sp = sp+hybrid(hyb,nums);
        end
        %
        rng(run);
        customers = customers(randperm(N));
        %
        taup = zeros(1,N); lambdap = zeros(1,N); mup = zeros(1,N);
        taub = zeros(1,N); lambdab = zeros(1,N); mub = zeros(1,N);
        for cs = 1:N
            if customers(cs) == 1
                taup(cs) = beta(h)/N; lambdap(cs) = beta(h)/N; mup(cs) = beta(h)/N;
                taub(cs) = beta(h)/N; lambdab(cs) = beta(h)/N; mub(cs) = beta(h)/N;
            elseif customers(cs) == 2
                taup(cs) = beta(h)/N; lambdap(cs) = 2*beta(h)/N; mup(cs) = beta(h)/N;
                taub(cs) = beta(h)/N; lambdab(cs) = 2*beta(h)/N; mub(cs) = beta(h)/N;
            elseif customers(cs) == 3
                taup(cs) = kappa/N; lambdap(cs) = 2*kappa/N; mup(cs) = 0;
                taub(cs) = kappa/N; lambdab(cs) = 2*kappa/N; mub(cs) = 0;
            elseif customers(cs) == 4
                taup(cs) = beta(h)/N; mup(cs) = beta(h)/N;
                lambdap(cs) = (N + nr_neighbors(cs))*beta(h)/(N^2);
                taub(cs) = beta(h)/N; mub(cs) = beta(h)/N;
                lambdab(cs) = (N + nr_neighbors(cs))*beta(h)/(N^2);
            end
        end
        rhop = 1./(2.*(taup+e));
        a_HP = rhop./(1 + lambdap.*rhop.*sigma_ij.*(N-1));
        b_HP = rhop./(1 + lambdap.*rhop.*(N-1));
        rhob = 1./(2.*(taub+e));
        a_HB = rhob./(1 - lambdab.*rhob);
        b_HB = rhob./(1 + lambdab.*rhob.*(N-1));
        
        for h = 1:Hours
            k=0;
            if time_zone(1,h) == 1
                k = 1;
            end
            if hyb == 1
               consumption_H1(run,1,:,h) = b_HP'.*(g_bar(k)-omega_bar.*mup)' + a_HP'.*(g(k,:) - g_bar(k))';
               total_consumption_H1(1,run,h) = sum(consumption_H1(run,1,:,h));
               price(h,:) = beta(h) * (total_consumption_H1(1,run,h) + omega)/N;
               price_container(run,13,h,:) = price(h,:);
               cost(h) = kappa * total_consumption_H1(1,run,h) ...
                   * total_consumption_H1(1,run,h)/N;
               revenue(h) = price(h,1)* total_consumption_H1(1,run,h);
               realized_profit_rate(run,13,h) = revenue(h)/ cost(h);
               cost_container_H1(run,1,h,1) = cost(h);
               revenue_container_H1(run,1,h) = revenue(h);
               for i = 1:N
                   u_i = -consumption_H1(run,1,i,h)*(beta(h)*(total_consumption_H1(1,run,h)...
                       +omega)/N)+g(k,i)*consumption_H1(run,1,i,h)-e(i)*consumption_H1(run,1,i,h)...
                       *consumption_H1(run,1,i,h);
                   utility_H1(run,1,i,h) = u_i;
                   welfare_H1(run,1,h) = welfare_H1(run,1,h) + u_i;
               end

               consumption_H1(run,2,:,h) = a_HB.*(g(k,:)-sum(g(k,:))/N) + b_HB.*(sum(g(k,:))/N-omega_bar.*mub);
               total_consumption_H1(2,run,h) = sum(consumption_H1(run,2,:,h));
               price(h,:) = beta(h) * (total_consumption_H1(2,run,h) + omega)/N;
               price_container(run,14,h,:) = price(h,:);
               cost(h) = kappa * total_consumption_H1(2,run,h) ...
                   * total_consumption_H1(2,run,h)/N;
               revenue(h) = price(h,1)* total_consumption_H1(2,run,h);
               realized_profit_rate(run,14,h) = revenue(h)/ cost(h);
               cost_container_H1(run,2,h,1) = cost(h);
               revenue_container_H1(run,2,h) = revenue(h);
               for i = 1:N
                   u_i = -consumption_H1(run,2,i,h)*(beta(h)*(total_consumption_H1(2,run,h)...
                       +omega)/N)+g(k,i)*consumption_H1(run,2,i,h)-e(i)*consumption_H1(run,2,i,h)...
                       *consumption_H1(run,2,i,h);
                   utility_H1(run,2,i,h) = u_i;
                   welfare_H1(run,2,h) = welfare_H1(run,2,h) + u_i;
               end
            elseif hyb == 2
               consumption_H2(run,1,:,h) = b_HP'.*(g_bar(k)-omega_bar.*mup)' + a_HP'.*(g(k,:) - g_bar(k))';
               total_consumption_H2(1,run,h) = sum(consumption_H2(run,1,:,h));
               price(h,:) = beta(h) * (total_consumption_H2(1,run,h) + omega)/N;
               price_container(run,16,h,:) = price(h,:);
               cost(h) = kappa * total_consumption_H2(1,run,h) ...
                   * total_consumption_H2(1,run,h)/N;
               revenue(h) = price(h,1)* total_consumption_H2(1,run,h);
               realized_profit_rate(run,16,h) = revenue(h)/ cost(h);
               cost_container_H2(run,1,h,1) = cost(h);
               revenue_container_H2(run,1,h) = revenue(h);
               for i = 1:N
                   u_i = -consumption_H2(run,1,i,h)*(beta(h)*(total_consumption_H2(1,run,h)...
                       +omega)/N)+g(k,i)*consumption_H2(run,1,i,h)-e(i)*consumption_H2(run,1,i,h)...
                       *consumption_H2(run,1,i,h);
                   utility_H2(run,1,i,h) = u_i;
                   welfare_H2(run,1,h) = welfare_H2(run,1,h) + u_i;
               end

               consumption_H2(run,2,:,h) = a_HB.*(g(k,:)-sum(g(k,:))/N) + b_HB.*(sum(g(k,:))/N-omega_bar.*mub);
               total_consumption_H2(2,run,h) = sum(consumption_H2(run,2,:,h));
               price(h,:) = beta(h) * (total_consumption_H2(2,run,h) + omega)/N;
               price_container(run,17,h,:) = price(h,:);
               cost(h) = kappa * total_consumption_H2(2,run,h) ...
                   * total_consumption_H2(2,run,h)/N;
               revenue(h) = price(h,1)* total_consumption_H2(2,run,h);
               realized_profit_rate(run,17,h) = revenue(h)/ cost(h);
               cost_container_H2(run,2,h,1) = cost(h);
               revenue_container_H2(run,2,h) = revenue(h);
               for i = 1:N
                   u_i = -consumption_H2(run,2,i,h)*(beta(h)*(total_consumption_H2(2,run,h)...
                       +omega)/N)+g(k,i)*consumption_H2(run,2,i,h)-e(i)*consumption_H2(run,2,i,h)...
                       *consumption_H2(run,2,i,h);
                   utility_H2(run,2,i,h) = u_i;
                   welfare_H2(run,2,h) = welfare_H2(run,2,h) + u_i;
               end
            elseif hyb == 3
               consumption_H3(run,1,:,h) = b_HP'.*(g_bar(k)-omega_bar.*mup)' + a_HP'.*(g(k,:) - g_bar(k))';
               total_consumption_H3(1,run,h) = sum(consumption_H3(run,1,:,h));
               price(h,:) = beta(h) * (total_consumption_H3(1,run,h) + omega)/N;
               price_container(run,19,h,:) = price(h,:);
               cost(h) = kappa * total_consumption_H3(1,run,h) ...
                   * total_consumption_H3(1,run,h)/N;
               revenue(h) = price(h,1)* total_consumption_H3(1,run,h);
               realized_profit_rate(run,19,h) = revenue(h)/ cost(h);
               cost_container_H3(run,1,h,1) = cost(h);
               revenue_container_H3(run,1,h) = revenue(h);
               for i = 1:N
                   u_i = -consumption_H3(run,1,i,h)*(beta(h)*(total_consumption_H3(1,run,h)...
                       +omega)/N)+g(k,i)*consumption_H3(run,1,i,h)-e(i)*consumption_H3(run,1,i,h)...
                       *consumption_H3(run,1,i,h);
                   utility_H3(run,1,i,h) = u_i;
                   welfare_H3(run,1,h) = welfare_H3(run,1,h) + u_i;
               end

               consumption_H3(run,2,:,h) = a_HB.*(g(k,:)-sum(g(k,:))/N) + b_HB.*(sum(g(k,:))/N-omega_bar.*mub);
               total_consumption_H3(2,run,h) = sum(consumption_H3(run,2,:,h));
               price(h,:) = beta(h) * (total_consumption_H3(2,run,h) + omega)/N;
               price_container(run,20,h,:) = price(h,:);
               cost(h) = kappa * total_consumption_H3(2,run,h) ...
                   * total_consumption_H3(2,run,h)/N;
               revenue(h) = price(h,1)* total_consumption_H3(2,run,h);
               realized_profit_rate(run,20,h) = revenue(h)/ cost(h);
               cost_container_H3(run,2,h,1) = cost(h);
               revenue_container_H3(run,2,h) = revenue(h);
               for i = 1:N
                   u_i = -consumption_H3(run,2,i,h)*(beta(h)*(total_consumption_H3(2,run,h)...
                       +omega)/N)+g(k,i)*consumption_H3(run,2,i,h)-e(i)*consumption_H3(run,2,i,h)...
                       *consumption_H3(run,2,i,h);
                   utility_H3(run,2,i,h) = u_i;
                   welfare_H3(run,2,h) = welfare_H3(run,2,h) + u_i;
               end
            elseif hyb == 4
               consumption_H4(run,1,:,h) = b_HP'.*(g_bar(k)-omega_bar.*mup)' + a_HP'.*(g(k,:) - g_bar(k))';
               total_consumption_H4(1,run,h) = sum(consumption_H4(run,1,:,h));
               price(h,:) = beta(h) * (total_consumption_H4(1,run,h) + omega)/N;
               price_container(run,22,h,:) = price(h,:);
               cost(h) = kappa * total_consumption_H4(1,run,h) ...
                   * total_consumption_H4(1,run,h)/N;
               revenue(h) = price(h,1)* total_consumption_H4(1,run,h);
               realized_profit_rate(run,22,h) = revenue(h)/ cost(h);
               cost_container_H4(run,1,h,1) = cost(h);
               revenue_container_H4(run,1,h) = revenue(h);
               for i = 1:N
                   u_i = -consumption_H4(run,1,i,h)*(beta(h)*(total_consumption_H4(1,run,h)...
                       +omega)/N)+g(k,i)*consumption_H4(run,1,i,h)-e(i)*consumption_H4(run,1,i,h)...
                       *consumption_H4(run,1,i,h);
                   utility_H4(run,1,i,h) = u_i;
                   welfare_H4(run,1,h) = welfare_H4(run,1,h) + u_i;
               end

               consumption_H4(run,2,:,h) = a_HB.*(g(k,:)-sum(g(k,:))/N) + b_HB.*(sum(g(k,:))/N-omega_bar.*mub);
               total_consumption_H4(2,run,h) = sum(consumption_H4(run,2,:,h));
               price(h,:) = beta(h) * (total_consumption_H4(2,run,h) + omega)/N;
               price_container(run,23,h,:) = price(h,:);
               cost(h) = kappa * total_consumption_H4(2,run,h) ...
                   * total_consumption_H4(2,run,h)/N;
               revenue(h) = price(h,1)* total_consumption_H4(2,run,h);
               realized_profit_rate(run,23,h) = revenue(h)/ cost(h);
               cost_container_H4(run,2,h,1) = cost(h);
               revenue_container_H4(run,2,h) = revenue(h);
               for i = 1:N
                   u_i = -consumption_H4(run,2,i,h)*(beta(h)*(total_consumption_H4(2,run,h)...
                       +omega)/N)+g(k,i)*consumption_H4(run,2,i,h)-e(i)*consumption_H4(run,2,i,h)...
                       *consumption_H4(run,2,i,h);
                   utility_H4(run,2,i,h) = u_i;
                   welfare_H4(run,2,h) = welfare_H4(run,2,h) + u_i;
               end
            elseif hyb == 5
               consumption_H5(run,1,:,h) = b_HP'.*(g_bar(k)-omega_bar.*mup)' + a_HP'.*(g(k,:) - g_bar(k))';
               total_consumption_H5(1,run,h) = sum(consumption_H5(run,1,:,h));
               price(h,:) = beta(h) * (total_consumption_H5(1,run,h) + omega)/N;
               price_container(run,25,h,:) = price(h,:);
               cost(h) = kappa * total_consumption_H5(1,run,h) ...
                   * total_consumption_H5(1,run,h)/N;
               revenue(h) = price(h,1)* total_consumption_H5(1,run,h);
               realized_profit_rate(run,25,h) = revenue(h)/ cost(h);
               cost_container_H5(run,1,h,1) = cost(h);
               revenue_container_H5(run,1,h) = revenue(h);
               for i = 1:N
                   u_i = -consumption_H5(run,1,i,h)*(beta(h)*(total_consumption_H5(1,run,h)...
                       +omega)/N)+g(k,i)*consumption_H5(run,1,i,h)-e(i)*consumption_H5(run,1,i,h)...
                       *consumption_H5(run,1,i,h);
                   utility_H5(run,1,i,h) = u_i;
                   welfare_H5(run,1,h) = welfare_H5(run,1,h) + u_i;
               end

               consumption_H5(run,2,:,h) = a_HB.*(g(k,:)-sum(g(k,:))/N) + b_HB.*(sum(g(k,:))/N-omega_bar.*mub);
               total_consumption_H5(2,run,h) = sum(consumption_H5(run,2,:,h));
               price(h,:) = beta(h) * (total_consumption_H5(2,run,h) + omega)/N;
               price_container(run,26,h,:) = price(h,:);
               cost(h) = kappa * total_consumption_H5(2,run,h) ...
                   * total_consumption_H5(2,run,h)/N;
               revenue(h) = price(h,1)* total_consumption_H5(2,run,h);
               realized_profit_rate(run,26,h) = revenue(h)/ cost(h);
               cost_container_H5(run,2,h,1) = cost(h);
               revenue_container_H5(run,2,h) = revenue(h);
               for i = 1:N
                   u_i = -consumption_H5(run,2,i,h)*(beta(h)*(total_consumption_H5(2,run,h)...
                       +omega)/N)+g(k,i)*consumption_H5(run,2,i,h)-e(i)*consumption_H5(run,2,i,h)...
                       *consumption_H5(run,2,i,h);
                   utility_H5(run,2,i,h) = u_i;
                   welfare_H5(run,2,h) = welfare_H5(run,2,h) + u_i;
               end
            elseif hyb == 6
               consumption_H6(run,1,:,h) = b_HP'.*(g_bar(k)-omega_bar.*mup)' + a_HP'.*(g(k,:) - g_bar(k))';
               total_consumption_H6(1,run,h) = sum(consumption_H6(run,1,:,h));
               price(h,:) = beta(h) * (total_consumption_H6(1,run,h) + omega)/N;
               price_container(run,28,h,:) = price(h,:);
               cost(h) = kappa * total_consumption_H6(1,run,h) ...
                   * total_consumption_H6(1,run,h)/N;
               revenue(h) = price(h,1)* total_consumption_H6(1,run,h);
               realized_profit_rate(run,28,h) = revenue(h)/ cost(h);
               cost_container_H6(run,1,h,1) = cost(h);
               revenue_container_H6(run,1,h) = revenue(h);
               for i = 1:N
                   u_i = -consumption_H6(run,1,i,h)*(beta(h)*(total_consumption_H6(1,run,h)...
                       +omega)/N)+g(k,i)*consumption_H6(run,1,i,h)-e(i)*consumption_H6(run,1,i,h)...
                       *consumption_H6(run,1,i,h);
                   utility_H6(run,1,i,h) = u_i;
                   welfare_H6(run,1,h) = welfare_H6(run,1,h) + u_i;
               end

               consumption_H6(run,2,:,h) = a_HB.*(g(k,:)-sum(g(k,:))/N) + b_HB.*(sum(g(k,:))/N-omega_bar.*mub);
               total_consumption_H6(2,run,h) = sum(consumption_H6(run,2,:,h));
               price(h,:) = beta(h) * (total_consumption_H6(2,run,h) + omega)/N;
               price_container(run,29,h,:) = price(h,:);
               cost(h) = kappa * total_consumption_H6(2,run,h) ...
                   * total_consumption_H6(2,run,h)/N;
               revenue(h) = price(h,1)* total_consumption_H6(2,run,h);
               realized_profit_rate(run,29,h) = revenue(h)/ cost(h);
               cost_container_H6(run,2,h,1) = cost(h);
               revenue_container_H6(run,2,h) = revenue(h);
               for i = 1:N
                   u_i = -consumption_H6(run,2,i,h)*(beta(h)*(total_consumption_H6(2,run,h)...
                       +omega)/N)+g(k,i)*consumption_H6(run,2,i,h)-e(i)*consumption_H6(run,2,i,h)...
                       *consumption_H6(run,2,i,h);
                   utility_H6(run,2,i,h) = u_i;
                   welfare_H6(run,2,h) = welfare_H6(run,2,h) + u_i;
               end
            elseif hyb == 7
               consumption_H7(run,1,:,h) = b_HP'.*(g_bar(k)-omega_bar.*mup)' + a_HP'.*(g(k,:) - g_bar(k))';
               total_consumption_H7(1,run,h) = sum(consumption_H7(run,1,:,h));
               price(h,:) = beta(h) * (total_consumption_H7(1,run,h) + omega)/N;
               price_container(run,31,h,:) = price(h,:);
               cost(h) = kappa * total_consumption_H7(1,run,h) ...
                   * total_consumption_H7(1,run,h)/N;
               revenue(h) = price(h,1)* total_consumption_H7(1,run,h);
               realized_profit_rate(run,31,h) = revenue(h)/ cost(h);
               cost_container_H7(run,1,h,1) = cost(h);
               revenue_container_H7(run,1,h) = revenue(h);
               for i = 1:N
                   u_i = -consumption_H7(run,1,i,h)*(beta(h)*(total_consumption_H7(1,run,h)...
                       +omega)/N)+g(k,i)*consumption_H7(run,1,i,h)-e(i)*consumption_H7(run,1,i,h)...
                       *consumption_H7(run,1,i,h);
                   utility_H7(run,1,i,h) = u_i;
                   welfare_H7(run,1,h) = welfare_H7(run,1,h) + u_i;
               end

               consumption_H7(run,2,:,h) = a_HB.*(g(k,:)-sum(g(k,:))/N) + b_HB.*(sum(g(k,:))/N-omega_bar.*mub);
               total_consumption_H7(2,run,h) = sum(consumption_H7(run,2,:,h));
               price(h,:) = beta(h) * (total_consumption_H7(2,run,h) + omega)/N;
               price_container(run,32,h,:) = price(h,:);
               cost(h) = kappa * total_consumption_H7(2,run,h) ...
                   * total_consumption_H7(2,run,h)/N;
               revenue(h) = price(h,1)* total_consumption_H7(2,run,h);
               realized_profit_rate(run,32,h) = revenue(h)/ cost(h);
               cost_container_H7(run,2,h,1) = cost(h);
               revenue_container_H7(run,2,h) = revenue(h);
               for i = 1:N
                   u_i = -consumption_H7(run,2,i,h)*(beta(h)*(total_consumption_H7(2,run,h)...
                       +omega)/N)+g(k,i)*consumption_H7(run,2,i,h)-e(i)*consumption_H7(run,2,i,h)...
                       *consumption_H7(run,2,i,h);
                   utility_H7(run,2,i,h) = u_i;
                   welfare_H7(run,2,h) = welfare_H7(run,2,h) + u_i;
               end
            elseif hyb == 8
               consumption_H8(run,1,:,h) = b_HP'.*(g_bar(k)-omega_bar.*mup)' + a_HP'.*(g(k,:) - g_bar(k))';
               total_consumption_H8(1,run,h) = sum(consumption_H8(run,1,:,h));
               price(h,:) = beta(h) * (total_consumption_H8(1,run,h) + omega)/N;
               price_container(run,34,h,:) = price(h,:);
               cost(h) = kappa * total_consumption_H8(1,run,h) ...
                   * total_consumption_H8(1,run,h)/N;
               revenue(h) = price(h,1)* total_consumption_H8(1,run,h);
               realized_profit_rate(run,34,h) = revenue(h)/ cost(h);
               cost_container_H8(run,1,h,1) = cost(h);
               revenue_container_H8(run,1,h) = revenue(h);
               for i = 1:N
                   u_i = -consumption_H8(run,1,i,h)*(beta(h)*(total_consumption_H8(1,run,h)...
                       +omega)/N)+g(k,i)*consumption_H8(run,1,i,h)-e(i)*consumption_H8(run,1,i,h)...
                       *consumption_H8(run,1,i,h);
                   utility_H8(run,1,i,h) = u_i;
                   welfare_H8(run,1,h) = welfare_H8(run,1,h) + u_i;
               end

               consumption_H8(run,2,:,h) = a_HB.*(g(k,:)-sum(g(k,:))/N) + b_HB.*(sum(g(k,:))/N-omega_bar.*mub);
               total_consumption_H8(2,run,h) = sum(consumption_H8(run,2,:,h));
               price(h,:) = beta(h) * (total_consumption_H8(2,run,h) + omega)/N;
               price_container(run,35,h,:) = price(h,:);
               cost(h) = kappa * total_consumption_H8(2,run,h) ...
                   * total_consumption_H8(2,run,h)/N;
               revenue(h) = price(h,1)* total_consumption_H8(2,run,h);
               realized_profit_rate(run,35,h) = revenue(h)/ cost(h);
               cost_container_H8(run,2,h,1) = cost(h);
               revenue_container_H8(run,2,h) = revenue(h);
               for i = 1:N
                   u_i = -consumption_H8(run,2,i,h)*(beta(h)*(total_consumption_H8(2,run,h)...
                       +omega)/N)+g(k,i)*consumption_H8(run,2,i,h)-e(i)*consumption_H8(run,2,i,h)...
                       *consumption_H8(run,2,i,h);
                   utility_H8(run,2,i,h) = u_i;
                   welfare_H8(run,2,h) = welfare_H8(run,2,h) + u_i;
               end
            elseif hyb == 9
               consumption_H9(run,1,:,h) = b_HP'.*(g_bar(k)-omega_bar.*mup)' + a_HP'.*(g(k,:) - g_bar(k))';
               total_consumption_H9(1,run,h) = sum(consumption_H9(run,1,:,h));
               price(h,:) = beta(h) * (total_consumption_H9(1,run,h) + omega)/N;
               price_container(run,37,h,:) = price(h,:);
               cost(h) = kappa * total_consumption_H9(1,run,h) ...
                   * total_consumption_H9(1,run,h)/N;
               revenue(h) = price(h,1)* total_consumption_H9(1,run,h);
               realized_profit_rate(run,37,h) = revenue(h)/ cost(h);
               cost_container_H9(run,1,h,1) = cost(h);
               revenue_container_H9(run,1,h) = revenue(h);
               for i = 1:N
                   u_i = -consumption_H9(run,1,i,h)*(beta(h)*(total_consumption_H9(1,run,h)...
                       +omega)/N)+g(k,i)*consumption_H9(run,1,i,h)-e(i)*consumption_H9(run,1,i,h)...
                       *consumption_H9(run,1,i,h);
                   utility_H9(run,1,i,h) = u_i;
                   welfare_H9(run,1,h) = welfare_H9(run,1,h) + u_i;
               end

               consumption_H9(run,2,:,h) = a_HB.*(g(k,:)-sum(g(k,:))/N) + b_HB.*(sum(g(k,:))/N-omega_bar.*mub);
               total_consumption_H9(2,run,h) = sum(consumption_H9(run,2,:,h));
               price(h,:) = beta(h) * (total_consumption_H9(2,run,h) + omega)/N;
               price_container(run,38,h,:) = price(h,:);
               cost(h) = kappa * total_consumption_H9(2,run,h) ...
                   * total_consumption_H9(2,run,h)/N;
               revenue(h) = price(h,1)* total_consumption_H9(2,run,h);
               realized_profit_rate(run,38,h) = revenue(h)/ cost(h);
               cost_container_H9(run,2,h,1) = cost(h);
               revenue_container_H9(run,2,h) = revenue(h);
               for i = 1:N
                   u_i = -consumption_H9(run,2,i,h)*(beta(h)*(total_consumption_H9(2,run,h)...
                       +omega)/N)+g(k,i)*consumption_H9(run,2,i,h)-e(i)*consumption_H9(run,2,i,h)...
                       *consumption_H9(run,2,i,h);
                   utility_H9(run,2,i,h) = u_i;
                   welfare_H9(run,2,h) = welfare_H9(run,2,h) + u_i;
               end
            elseif hyb == 10
               consumption_H10(run,1,:,h) = b_HP'.*(g_bar(k)-omega_bar.*mup)' + a_HP'.*(g(k,:) - g_bar(k))';
               total_consumption_H10(1,run,h) = sum(consumption_H10(run,1,:,h));
               price(h,:) = beta(h) * (total_consumption_H10(1,run,h) + omega)/N;
               price_container(run,40,h,:) = price(h,:);
               cost(h) = kappa * total_consumption_H10(1,run,h) ...
                   * total_consumption_H10(1,run,h)/N;
               revenue(h) = price(h,1)* total_consumption_H10(1,run,h);
               realized_profit_rate(run,40,h) = revenue(h)/ cost(h);
               cost_container_H10(run,1,h,1) = cost(h);
               revenue_container_H10(run,1,h) = revenue(h);
               for i = 1:N
                   u_i = -consumption_H10(run,1,i,h)*(beta(h)*(total_consumption_H10(1,run,h)...
                       +omega)/N)+g(k,i)*consumption_H10(run,1,i,h)-e(i)*consumption_H10(run,1,i,h)...
                       *consumption_H10(run,1,i,h);
                   utility_H10(run,1,i,h) = u_i;
                   welfare_H10(run,1,h) = welfare_H10(run,1,h) + u_i;
               end

               consumption_H10(run,2,:,h) = a_HB.*(g(k,:)-sum(g(k,:))/N) + b_HB.*(sum(g(k,:))/N-omega_bar.*mub);
               total_consumption_H10(2,run,h) = sum(consumption_H10(run,2,:,h));
               price(h,:) = beta(h) * (total_consumption_H10(2,run,h) + omega)/N;
               price_container(run,41,h,:) = price(h,:);
               cost(h) = kappa * total_consumption_H10(2,run,h) ...
                   * total_consumption_H10(2,run,h)/N;
               revenue(h) = price(h,1)* total_consumption_H10(2,run,h);
               realized_profit_rate(run,41,h) = revenue(h)/ cost(h);
               cost_container_H10(run,2,h,1) = cost(h);
               revenue_container_H10(run,2,h) = revenue(h);
               for i = 1:N
                   u_i = -consumption_H10(run,2,i,h)*(beta(h)*(total_consumption_H10(2,run,h)...
                       +omega)/N)+g(k,i)*consumption_H10(run,2,i,h)-e(i)*consumption_H10(run,2,i,h)...
                       *consumption_H10(run,2,i,h);
                   utility_H10(run,2,i,h) = u_i;
                   welfare_H10(run,2,h) = welfare_H10(run,2,h) + u_i;
               end
            end
        end 
    end
end
toc;

%% Plots
ccc = colormap(hsv(8));

average_consumption_S = reshape(mean(total_consumption_S,2),[3,6]);
average_consumption_U = reshape(mean(total_consumption_U,2),[3,6]);
average_consumption_W = reshape(mean(total_consumption_W,2),[3,6]);
average_consumption_C = reshape(mean(total_consumption_C,2),[3,6]);

fig = figure;
hold on
plot(1:1:Hours,average_consumption_S(1,:)/N,'s-','Color',ccc(1,:),'LineWidth',2 )
plot(1:1:Hours,average_consumption_U(1,:)/N,'s-','Color',ccc(3,:),'LineWidth',2 )
plot(1:1:Hours,average_consumption_W(1,:)/N,'s-','Color',ccc(5,:),'LineWidth',2 )
plot(1:1:Hours,average_consumption_C(1,:)/N,'s-','Color',ccc(7,:),'LineWidth',2 )
grid on
title('Private Information')
ylabel('Average Consumption per Capita')
xlabel('Hours')
lll = legend('$\Gamma = $S, $\Omega = $P','$\Gamma = $U, $\Omega = $P',...
    '$\Gamma = $W, $\Omega = $P','$\Gamma = $CF, $\Omega = $P');
set(lll,'Interpreter','Latex');
set(lll,'FontSize',14);
hold off

fig = figure;
hold on
plot(1:1:Hours,average_consumption_S(2,:)/N,'s-','Color',ccc(1,:),'LineWidth',2 )
plot(1:1:Hours,average_consumption_U(2,:)/N,'s-','Color',ccc(3,:),'LineWidth',2 )
plot(1:1:Hours,average_consumption_W(2,:)/N,'s-','Color',ccc(5,:),'LineWidth',2 )
plot(1:1:Hours,average_consumption_C(2,:)/N,'s-','Color',ccc(7,:),'LineWidth',2 )
grid on
title('Broadcast')
ylabel('Average Consumption per Capita')
xlabel('Hours')
lll = legend('$\Gamma = $S, $\Omega = $B','$\Gamma = $U, $\Omega = $B',...
    '$\Gamma = $W, $\Omega = $B','$\Gamma = $CF, $\Omega = $B');
set(lll,'Interpreter','Latex');
set(lll,'FontSize',14);
hold off

average_consumption_H1 = reshape(mean(total_consumption_H1,2),[3,6]);
average_consumption_H2 = reshape(mean(total_consumption_H2,2),[3,6]);
average_consumption_H3 = reshape(mean(total_consumption_H3,2),[3,6]);

fig = figure;
hold on
plot(1:1:Hours,average_consumption_H1(1,:)/N,'s-','Color',ccc(2,:),'LineWidth',2 )
plot(1:1:Hours,average_consumption_H2(1,:)/N,'s-','Color',ccc(4,:),'LineWidth',2 )
plot(1:1:Hours,average_consumption_H3(1,:)/N,'s-','Color',ccc(6,:),'LineWidth',2 )
grid on
title('Private Information - Hybrid Models with S')
ylabel('Average Consumption per Capita')
xlabel('Hours')
lll = legend('$\Gamma = $H1, $\Omega = $P','$\Gamma = $H2, $\Omega = $P',...
    '$\Gamma = $H3, $\Omega = $P');
set(lll,'Interpreter','Latex');
set(lll,'FontSize',14);
hold off

fig = figure;
hold on
plot(1:1:Hours,average_consumption_H1(2,:)/N,'s-','Color',ccc(2,:),'LineWidth',2 )
plot(1:1:Hours,average_consumption_H2(2,:)/N,'s-','Color',ccc(4,:),'LineWidth',2 )
plot(1:1:Hours,average_consumption_H3(2,:)/N,'s-','Color',ccc(6,:),'LineWidth',2 )
grid on
title('Broadcast - Hybrid Models with S')
ylabel('Average Consumption per Capita')
xlabel('Hours')
lll = legend('$\Gamma = $H1, $\Omega = $B','$\Gamma = $H2, $\Omega = $B',...
    '$\Gamma = $H3, $\Omega = $B');
set(lll,'Interpreter','Latex');
set(lll,'FontSize',14);
hold off

average_consumption_H4 = reshape(mean(total_consumption_H4,2),[3,6]);
average_consumption_H5 = reshape(mean(total_consumption_H5,2),[3,6]);
average_consumption_H6 = reshape(mean(total_consumption_H6,2),[3,6]);

fig = figure;
hold on
plot(1:1:Hours,average_consumption_H4(1,:)/N,'s-','Color',ccc(2,:),'LineWidth',2 )
plot(1:1:Hours,average_consumption_H5(1,:)/N,'s-','Color',ccc(4,:),'LineWidth',2 )
plot(1:1:Hours,average_consumption_H6(1,:)/N,'s-','Color',ccc(6,:),'LineWidth',2 )
grid on
title('Private Information - Hybrid Model with CF')
ylabel('Average Consumption per Capita')
xlabel('Hours')
lll = legend('$\Gamma = $H4, $\Omega = $P','$\Gamma = $H5, $\Omega = $P',...
    '$\Gamma = $H6, $\Omega = $P');
set(lll,'Interpreter','Latex');
set(lll,'FontSize',14);
hold off

fig = figure;
hold on
plot(1:1:Hours,average_consumption_H4(2,:)/N,'s-','Color',ccc(2,:),'LineWidth',2 )
plot(1:1:Hours,average_consumption_H5(2,:)/N,'s-','Color',ccc(4,:),'LineWidth',2 )
plot(1:1:Hours,average_consumption_H6(2,:)/N,'s-','Color',ccc(6,:),'LineWidth',2 )
grid on
title('Broadcast - Hybrid Models with CF')
ylabel('Average Consumption per Capita')
xlabel('Hours')
lll = legend('$\Gamma = $H4, $\Omega = $B','$\Gamma = $H5, $\Omega = $B',...
    '$\Gamma = $H6, $\Omega = $B');
set(lll,'Interpreter','Latex');
set(lll,'FontSize',14);
hold off

average_consumption_H7 = reshape(mean(total_consumption_H7,2),[3,6]);
fig = figure;
hold on
plot(1:1:Hours,average_consumption_H7(1,:)/N,'s-','Color',ccc(2,:),'LineWidth',2 )
plot(1:1:Hours,average_consumption_H7(2,:)/N,'s-','Color',ccc(8,:),'LineWidth',2 )
grid on
title('Hybrid Model - 7')
ylabel('Average Consumption per Capita')
xlabel('Hours')
lll = legend('$\Gamma = $H7, $\Omega = $P','$\Gamma = $H7, $\Omega = $B');
set(lll,'Interpreter','Latex');
set(lll,'FontSize',14);
hold off

average_consumption_H8 = reshape(mean(total_consumption_H8,2),[3,6]);
fig = figure;
hold on
plot(1:1:Hours,average_consumption_H8(1,:)/N,'s-','Color',ccc(2,:),'LineWidth',2 )
plot(1:1:Hours,average_consumption_H8(2,:)/N,'s-','Color',ccc(8,:),'LineWidth',2 )
grid on
title('Hybrid Model - 8')
ylabel('Average Consumption per Capita')
xlabel('Hours')
lll = legend('$\Gamma = $H8, $\Omega = $P','$\Gamma = $H8, $\Omega = $B');
set(lll,'Interpreter','Latex');
set(lll,'FontSize',14);
hold off

average_consumption_H9 = reshape(mean(total_consumption_H9,2),[3,6]);
fig = figure;
hold on
plot(1:1:Hours,average_consumption_H9(1,:)/N,'s-','Color',ccc(2,:),'LineWidth',2 )
plot(1:1:Hours,average_consumption_H9(2,:)/N,'s-','Color',ccc(8,:),'LineWidth',2 )
grid on
title('Hybrid Model - 9')
ylabel('Average Consumption per Capita')
xlabel('Hours')
lll = legend('$\Gamma = $H9, $\Omega = $P','$\Gamma = $H9, $\Omega = $B');
set(lll,'Interpreter','Latex');
set(lll,'FontSize',14);
hold off

average_consumption_H10 = reshape(mean(total_consumption_H10,2),[3,6]);
fig = figure;
hold on
plot(1:1:Hours,average_consumption_H10(1,:)/N,'s-','Color',ccc(2,:),'LineWidth',2 )
plot(1:1:Hours,average_consumption_H10(2,:)/N,'s-','Color',ccc(8,:),'LineWidth',2 )
grid on
title('Hybrid Model - 10')
ylabel('Average Consumption per Capita')
xlabel('Hours')
lll = legend('$\Gamma = $H10, $\Omega = $P','$\Gamma = $H10, $\Omega = $B');
set(lll,'Interpreter','Latex');
set(lll,'FontSize',14);
hold off

%% Save Data
% Average Consumptions - needs to be divided by N
save './Results_Graph/sim1_avcons_S' average_consumption_S
save './Results_Graph/sim1_avcons_U' average_consumption_U
save './Results_Graph/sim1_avcons_W' average_consumption_W
save './Results_Graph/sim1_avcons_C' average_consumption_C
save './Results_Graph/sim1_avcons_H1' average_consumption_H1
save './Results_Graph/sim1_avcons_H2' average_consumption_H2
save './Results_Graph/sim1_avcons_H3' average_consumption_H3
save './Results_Graph/sim1_avcons_H4' average_consumption_H4
save './Results_Graph/sim1_avcons_H5' average_consumption_H5
save './Results_Graph/sim1_avcons_H6' average_consumption_H6
save './Results_Graph/sim1_avcons_H7' average_consumption_H7
save './Results_Graph/sim1_avcons_H8' average_consumption_H8
save './Results_Graph/sim1_avcons_H9' average_consumption_H9
save './Results_Graph/sim1_avcons_H10' average_consumption_H10
% Utility Containers
save './Results_Graph/sim1_utility_S' utility_S
save './Results_Graph/sim1_utility_U' utility_U
save './Results_Graph/sim1_utility_W' utility_W
save './Results_Graph/sim1_utility_C' utility_C
save './Results_Graph/sim1_utility_H1' utility_H1
save './Results_Graph/sim1_utility_H2' utility_H2
save './Results_Graph/sim1_utility_H3' utility_H3
save './Results_Graph/sim1_utility_H4' utility_H4
save './Results_Graph/sim1_utility_H5' utility_H5
save './Results_Graph/sim1_utility_H6' utility_H6
save './Results_Graph/sim1_utility_H7' utility_H7
save './Results_Graph/sim1_utility_H8' utility_H8
save './Results_Graph/sim1_utility_H9' utility_H9
save './Results_Graph/sim1_utility_H10' utility_H10
% Price Container
save './Results_Graph/sim1_price_container' price_container
% Revenue Containers
save './Results_Graph/sim1_revenue_container_S' revenue_container_S
save './Results_Graph/sim1_revenue_container_U' revenue_container_U
save './Results_Graph/sim1_revenue_container_W' revenue_container_W
save './Results_Graph/sim1_revenue_container_C' revenue_container_C
save './Results_Graph/sim1_revenue_container_H1' revenue_container_H1
save './Results_Graph/sim1_revenue_container_H2' revenue_container_H2
save './Results_Graph/sim1_revenue_container_H3' revenue_container_H3
save './Results_Graph/sim1_revenue_container_H4' revenue_container_H4
save './Results_Graph/sim1_revenue_container_H5' revenue_container_H5
save './Results_Graph/sim1_revenue_container_H6' revenue_container_H6
save './Results_Graph/sim1_revenue_container_H7' revenue_container_H7
save './Results_Graph/sim1_revenue_container_H8' revenue_container_H8
save './Results_Graph/sim1_revenue_container_H7' revenue_container_H9
save './Results_Graph/sim1_revenue_container_H8' revenue_container_H10
% Cost Containers
save './Results_Graph/sim1_cost_container_S' cost_container_S
save './Results_Graph/sim1_cost_container_U' cost_container_U
save './Results_Graph/sim1_cost_container_W' cost_container_W
save './Results_Graph/sim1_cost_container_C' cost_container_C
save './Results_Graph/sim1_cost_container_H1' cost_container_H1
save './Results_Graph/sim1_cost_container_H2' cost_container_H2
save './Results_Graph/sim1_cost_container_H3' cost_container_H3
save './Results_Graph/sim1_cost_container_H4' cost_container_H4
save './Results_Graph/sim1_cost_container_H5' cost_container_H5
save './Results_Graph/sim1_cost_container_H6' cost_container_H6
save './Results_Graph/sim1_cost_container_H7' cost_container_H7
save './Results_Graph/sim1_cost_container_H8' cost_container_H8
save './Results_Graph/sim1_cost_container_H9' cost_container_H9
save './Results_Graph/sim1_cost_container_H10' cost_container_H10

%% Simulation - 24 Hours
tic;

% Data Cleaning
clear beta_container

% Initilization
Hours = 24;

% Total Consumption Matricies
total_consumption_24 = zeros(8,run_time,Hours);

% Estimated and Realized Value Matricies
estimated_total_consumption_container = zeros(run_time,Hours);
estimated_total_revenue_container = zeros(run_time,Hours);
estimated_total_cost_container = zeros(run_time,Hours);
estimated_profit_container = zeros(run_time,1);
estimated_welfare_container = zeros(run_time,Hours);

realized_profit_rate_24 = zeros(run_time,8,Hours);
price_container_24 = zeros(run_time,8,Hours,N);

% Cost Containers
cost_container_24 = zeros(run_time,8,Hours,N);

% Revenue Containers
revenue_container_24 = zeros(run_time,8,Hours);

% Welfare Containers
welfare_24 = zeros(run_time,8,Hours);

% Consumption
consumption_24 = zeros(run_time,8,N,Hours);

% Utility Container
utility_24 = zeros(run_time,8,N,Hours);

% Other Matricies
estimated_revenue = zeros(run_time,Hours);
estimated_cost = zeros(run_time,Hours);
estimated_consumption = zeros(run_time,N,Hours);
estimated_total_consumption = zeros(run_time,Hours);
expected_total_consumption_squared = zeros(run_time,Hours);

% Other Parameters and Configurations
realized_omega = zeros(run_time,1);
seed= rng;
rng('default')
sigma_pop = 1;
sigma_ij = 0;

for run = 1:run_time
    current_sigma_ij = sigma_ij;
    time_zone = ones(1,Hours);
    
    graph = generate_graph(N,run);
    incidence_matrix = graph;
    nr_nodes = N;
    nr_neighbors = incidence_matrix*ones(nr_nodes,1);
    
    beta = zeros(Hours,1);
    
    % Setup
    g_bar = mean(cons_pref_rand);
    g = zeros(1,N);
    e = zeros(1,N);
    sigma_ii = 4;
    Sigma = sigma_ii * eye(N,N) + current_sigma_ij * (ones(N,N)-eye(N,N));
    
    % Inference Matrix
    IM = zeros(N,N);
    for i = 1:N
        for j = 1:N
            if i ~= j
                IM(i,j) = Sigma(i,j)/Sigma(i,i);
            end
        end
    end
    
    omega_bar = 0;
    omega_variance = 2;
    profit_rate = [1.0; 1.2; 1.1; 1.3]; %Fixed profit rate
    profit_rate = repelem(profit_rate,6);
    kappa = 1;
    beta(:) = profit_rate .* kappa .* ones(Hours,1);
    beta_container(:,run) = beta(:);
    sigma_ij = current_sigma_ij/sigma_ii;
    
    % Estimation of consumption, cost etc.
    for h = 1:Hours
       k = 0;
       if time_zone(1,h) == 1
           k = 1;
       end
       a = ones(N,1)/((N+1)*beta(h) + 2*mean(cons_pref_rand));
       rho = 1/(2*(beta(h) + mean(cons_pref_rand)));
       b = ones(N,1)*rho/((N-1)*beta(h)*rho*sigma_ij+1);

       estimated_consumption(run,:,h) = a*(g_bar(k)-omega_bar*beta(h));
       estimated_total_consumption(run,h) = sum(estimated_consumption(run,:,h));
       expected_total_consumption_squared(run,h) = estimated_total_consumption(run,h)...
           *estimated_total_consumption(run,h) + sigma_ii*sigma_ii* b'*b ;
       estimated_revenue(run,h) = beta(h) * expected_total_consumption_squared(run,h) ...
           + omega_bar *estimated_total_consumption(run,h);
       estimated_cost(run,h) = kappa *expected_total_consumption_squared(run,h);
       expected_profit_rate(run,1) = estimated_revenue(run,h)/estimated_cost(run,h);
    end
    
    estimated_total_consumption_container(run,:) =  estimated_total_consumption(run,:);
    estimated_total_revenue_container(run,:) = estimated_revenue(run,:);
    estimated_total_cost_container(run,:) = estimated_cost(run,:);
    estimated_profit_container(run) = expected_profit_rate(run);
    total_estimated_cost = sum(estimated_cost);
    
    % Simulation
    rng(run)
    rand_indices = ceil(length(cons_pref_rand)*rand(1,N));
    g(1,:) = cons_pref_rand(rand_indices);
    e(1,:) = price_sensitivity_rand(rand_indices);
    
    omega = omega_bar + omega_variance * randn(1);
    realized_omega(run) = omega;
    cost = zeros(Hours,1);
    revenue = zeros(Hours,1);
    price = zeros(Hours,N);
    
    for hyb = 7:10
        for h = 1:Hours
            customers = zeros(1,N);
            sp = 1;
            for nums = 1:4
                if sp > 100
                    break
                end
                customers(sp:sp+hybrid(hyb,nums)-1) = nums;
                sp = sp+hybrid(hyb,nums);
            end
            %
            rng(run);
            customers = customers(randperm(N));
            %
            taup = zeros(1,N); lambdap = zeros(1,N); mup = zeros(1,N);
            taub = zeros(1,N); lambdab = zeros(1,N); mub = zeros(1,N);
            for cs = 1:N
                if customers(cs) == 1
                    taup(cs) = beta(h)/N; lambdap(cs) = beta(h)/N; mup(cs) = beta(h)/N;
                    taub(cs) = beta(h)/N; lambdab(cs) = beta(h)/N; mub(cs) = beta(h)/N;
                elseif customers(cs) == 2
                    taup(cs) = beta(h)/N; lambdap(cs) = 2*beta(h)/N; mup(cs) = beta(h)/N;
                    taub(cs) = beta(h)/N; lambdab(cs) = 2*beta(h)/N; mub(cs) = beta(h)/N;
                elseif customers(cs) == 3
                    taup(cs) = kappa/N; lambdap(cs) = 2*kappa/N; mup(cs) = 0;
                    taub(cs) = kappa/N; lambdab(cs) = 2*kappa/N; mub(cs) = 0;
                elseif customers(cs) == 4
                    taup(cs) = beta(h)/N; mup(cs) = beta(h)/N;
                    lambdap(cs) = (N + nr_neighbors(cs))*beta(h)/(N^2);
                    taub(cs) = beta(h)/N; mub(cs) = beta(h)/N;
                    lambdab(cs) = (N + nr_neighbors(cs))*beta(h)/(N^2);
                end
            end
            rhop = 1./(2.*(taup+e));
            a_HP = rhop./(1 + lambdap.*rhop.*sigma_ij.*(N-1));
            b_HP = rhop./(1 + lambdap.*rhop.*(N-1));
            rhob = 1./(2.*(taub+e));
            a_HB = rhob./(1 - lambdab.*rhob);
            b_HB = rhob./(1 + lambdab.*rhob.*(N-1));
        
            k=0;
            if time_zone(1,h) == 1
                k = 1;
            end
            if hyb == 7
               consumption_24(run,1,:,h) = b_HP'.*(g_bar(k)-omega_bar.*mup)' + a_HP'.*(g(k,:) - g_bar(k))';
               total_consumption_24(1,run,h) = sum(consumption_24(run,1,:,h));
               price(h,:) = beta(h) * (total_consumption_24(1,run,h) + omega)/N;
               price_container_24(run,1,h,:) = price(h,:);
               cost(h) = kappa * total_consumption_24(1,run,h) ...
                   * total_consumption_24(1,run,h)/N;
               revenue(h) = price(h,1)* total_consumption_24(1,run,h);
               realized_profit_rate_24(run,1,h) = revenue(h)/ cost(h);
               cost_container_24(run,1,h,1) = cost(h);
               revenue_container_24(run,1,h) = revenue(h);
               for i = 1:N
                   u_i = -consumption_24(run,1,i,h)*(beta(h)*(total_consumption_24(1,run,h)...
                       +omega)/N)+g(k,i)*consumption_24(run,1,i,h)-e(i)*consumption_24(run,1,i,h)...
                       *consumption_24(run,1,i,h);
                   utility_24(run,1,i,h) = u_i;
                   welfare_24(run,1,h) = welfare_24(run,1,h) + u_i;
               end

               consumption_24(run,2,:,h) = a_HB.*(g(k,:)-sum(g(k,:))/N) + b_HB.*(sum(g(k,:))/N-omega_bar.*mub);
               total_consumption_24(2,run,h) = sum(consumption_24(run,2,:,h));
               price(h,:) = beta(h) * (total_consumption_24(2,run,h) + omega)/N;
               price_container_24(run,2,h,:) = price(h,:);
               cost(h) = kappa * total_consumption_24(2,run,h) ...
                   * total_consumption_24(2,run,h)/N;
               revenue(h) = price(h,1)* total_consumption_24(2,run,h);
               realized_profit_rate_24(run,2,h) = revenue(h)/ cost(h);
               cost_container_24(run,2,h,1) = cost(h);
               revenue_container_24(run,2,h) = revenue(h);
               for i = 1:N
                   u_i = -consumption_24(run,2,i,h)*(beta(h)*(total_consumption_24(2,run,h)...
                       +omega)/N)+g(k,i)*consumption_24(run,2,i,h)-e(i)*consumption_24(run,2,i,h)...
                       *consumption_24(run,2,i,h);
                   utility_24(run,2,i,h) = u_i;
                   welfare_24(run,2,h) = welfare_24(run,2,h) + u_i;
               end
            elseif hyb == 8
               consumption_24(run,3,:,h) = b_HP'.*(g_bar(k)-omega_bar.*mup)' + a_HP'.*(g(k,:) - g_bar(k))';
               total_consumption_24(3,run,h) = sum(consumption_24(run,3,:,h));
               price(h,:) = beta(h) * (total_consumption_24(3,run,h) + omega)/N;
               price_container_24(run,3,h,:) = price(h,:);
               cost(h) = kappa * total_consumption_24(3,run,h) ...
                   * total_consumption_24(3,run,h)/N;
               revenue(h) = price(h,1)* total_consumption_24(3,run,h);
               realized_profit_rate_24(run,3,h) = revenue(h)/ cost(h);
               cost_container_24(run,3,h,1) = cost(h);
               revenue_container_24(run,3,h) = revenue(h);
               for i = 1:N
                   u_i = -consumption_24(run,3,i,h)*(beta(h)*(total_consumption_24(3,run,h)...
                       +omega)/N)+g(k,i)*consumption_24(run,3,i,h)-e(i)*consumption_24(run,3,i,h)...
                       *consumption_24(run,3,i,h);
                   utility_24(run,3,i,h) = u_i;
                   welfare_24(run,3,h) = welfare_24(run,3,h) + u_i;
               end

               consumption_24(run,4,:,h) = a_HB.*(g(k,:)-sum(g(k,:))/N) + b_HB.*(sum(g(k,:))/N-omega_bar.*mub);
               total_consumption_24(4,run,h) = sum(consumption_24(run,4,:,h));
               price(h,:) = beta(h) * (total_consumption_24(4,run,h) + omega)/N;
               price_container_24(run,4,h,:) = price(h,:);
               cost(h) = kappa * total_consumption_24(4,run,h) ...
                   * total_consumption_24(4,run,h)/N;
               revenue(h) = price(h,1)* total_consumption_24(4,run,h);
               realized_profit_rate_24(run,4,h) = revenue(h)/ cost(h);
               cost_container_24(run,4,h,1) = cost(h);
               revenue_container_24(run,4,h) = revenue(h);
               for i = 1:N
                   u_i = -consumption_24(run,4,i,h)*(beta(h)*(total_consumption_24(4,run,h)...
                       +omega)/N)+g(k,i)*consumption_24(run,4,i,h)-e(i)*consumption_24(run,4,i,h)...
                       *consumption_24(run,4,i,h);
                   utility_24(run,4,i,h) = u_i;
                   welfare_24(run,4,h) = welfare_24(run,4,h) + u_i;
               end
            elseif hyb == 9
               consumption_24(run,5,:,h) = b_HP'.*(g_bar(k)-omega_bar.*mup)' + a_HP'.*(g(k,:) - g_bar(k))';
               total_consumption_24(5,run,h) = sum(consumption_24(run,5,:,h));
               price(h,:) = beta(h) * (total_consumption_24(5,run,h) + omega)/N;
               price_container_24(run,5,h,:) = price(h,:);
               cost(h) = kappa * total_consumption_24(5,run,h) ...
                   * total_consumption_24(5,run,h)/N;
               revenue(h) = price(h,1)* total_consumption_24(5,run,h);
               realized_profit_rate_24(run,5,h) = revenue(h)/ cost(h);
               cost_container_24(run,5,h,1) = cost(h);
               revenue_container_24(run,5,h) = revenue(h);
               for i = 1:N
                   u_i = -consumption_24(run,5,i,h)*(beta(h)*(total_consumption_24(5,run,h)...
                       +omega)/N)+g(k,i)*consumption_24(run,5,i,h)-e(i)*consumption_24(run,5,i,h)...
                       *consumption_24(run,5,i,h);
                   utility_24(run,5,i,h) = u_i;
                   welfare_24(run,5,h) = welfare_24(run,5,h) + u_i;
               end

               consumption_24(run,6,:,h) = a_HB.*(g(k,:)-sum(g(k,:))/N) + b_HB.*(sum(g(k,:))/N-omega_bar.*mub);
               total_consumption_24(6,run,h) = sum(consumption_24(run,6,:,h));
               price(h,:) = beta(h) * (total_consumption_24(6,run,h) + omega)/N;
               price_container_24(run,6,h,:) = price(h,:);
               cost(h) = kappa * total_consumption_24(6,run,h) ...
                   * total_consumption_24(6,run,h)/N;
               revenue(h) = price(h,1)* total_consumption_24(6,run,h);
               realized_profit_rate_24(run,6,h) = revenue(h)/ cost(h);
               cost_container_24(run,6,h,1) = cost(h);
               revenue_container_24(run,6,h) = revenue(h);
               for i = 1:N
                   u_i = -consumption_24(run,6,i,h)*(beta(h)*(total_consumption_24(6,run,h)...
                       +omega)/N)+g(k,i)*consumption_24(run,6,i,h)-e(i)*consumption_24(run,6,i,h)...
                       *consumption_24(run,6,i,h);
                   utility_24(run,6,i,h) = u_i;
                   welfare_24(run,6,h) = welfare_24(run,6,h) + u_i;
               end
            elseif hyb == 10
               consumption_24(run,7,:,h) = b_HP'.*(g_bar(k)-omega_bar.*mup)' + a_HP'.*(g(k,:) - g_bar(k))';
               total_consumption_24(7,run,h) = sum(consumption_24(run,7,:,h));
               price(h,:) = beta(h) * (total_consumption_24(7,run,h) + omega)/N;
               price_container_24(run,7,h,:) = price(h,:);
               cost(h) = kappa * total_consumption_24(7,run,h) ...
                   * total_consumption_24(7,run,h)/N;
               revenue(h) = price(h,1)* total_consumption_24(7,run,h);
               realized_profit_rate_24(run,7,h) = revenue(h)/ cost(h);
               cost_container_24(run,7,h,1) = cost(h);
               revenue_container_24(run,7,h) = revenue(h);
               for i = 1:N
                   u_i = -consumption_24(run,7,i,h)*(beta(h)*(total_consumption_24(7,run,h)...
                       +omega)/N)+g(k,i)*consumption_24(run,7,i,h)-e(i)*consumption_24(run,7,i,h)...
                       *consumption_24(run,7,i,h);
                   utility_24(run,7,i,h) = u_i;
                   welfare_24(run,7,h) = welfare_24(run,7,h) + u_i;
               end

               consumption_24(run,8,:,h) = a_HB.*(g(k,:)-sum(g(k,:))/N) + b_HB.*(sum(g(k,:))/N-omega_bar.*mub);
               total_consumption_24(8,run,h) = sum(consumption_24(run,8,:,h));
               price(h,:) = beta(h) * (total_consumption_24(8,run,h) + omega)/N;
               price_container_24(run,8,h,:) = price(h,:);
               cost(h) = kappa * total_consumption_24(8,run,h) ...
                   * total_consumption_24(8,run,h)/N;
               revenue(h) = price(h,1)* total_consumption_24(8,run,h);
               realized_profit_rate_24(run,8,h) = revenue(h)/ cost(h);
               cost_container_24(run,8,h,1) = cost(h);
               revenue_container_24(run,8,h) = revenue(h);
               for i = 1:N
                   u_i = -consumption_24(run,8,i,h)*(beta(h)*(total_consumption_24(8,run,h)...
                       +omega)/N)+g(k,i)*consumption_24(run,8,i,h)-e(i)*consumption_24(run,8,i,h)...
                       *consumption_24(run,8,i,h);
                   utility_24(run,8,i,h) = u_i;
                   welfare_24(run,8,h) = welfare_24(run,8,h) + u_i;
               end
            end
        end 
    end
end

toc;

%% Save Data - 24 Hours
% Average Consumptions - needs to be divided by N
average_consumption_24 = reshape(mean(total_consumption_24,2),[8,24]);
save './Results_Graph/sim1_avcons_24' average_consumption_24
% Utility Containers
save './Results_Graph/sim1_utility_24' utility_24
% Price Container
save './Results_Graph/sim1_price_container_24' price_container_24
% Revenue Containers
save './Results_Graph/sim1_revenue_container_24' revenue_container_24
% Cost Containers
save './Results_Graph/sim1_cost_container_24' cost_container_24

%% Plots - 24 Hours
fig = figure;
hold on
plot(1:1:Hours,average_consumption_24(1,:)/N,'s-','Color',ccc(1,:),'LineWidth',2 )
plot(1:1:Hours,average_consumption_24(3,:)/N,'s-','Color',ccc(3,:),'LineWidth',2 )
plot(1:1:Hours,average_consumption_24(5,:)/N,'s-','Color',ccc(5,:),'LineWidth',2 )
plot(1:1:Hours,average_consumption_24(7,:)/N,'s-','Color',ccc(7,:),'LineWidth',2 )
grid on
title('Hybrid Models - 24 Hours')
ylabel('Average Consumption per Capita')
xlabel('Hours')
lll = legend('$\Gamma = $H7, $\Omega = $P','$\Gamma = $H8, $\Omega = $P',...
    '$\Gamma = $H9, $\Omega = $P','$\Gamma = $H10, $\Omega = $P');
set(lll,'Interpreter','Latex');
set(lll,'FontSize',14);
hold off

fig = figure;
hold on
plot(1:1:Hours,average_consumption_24(2,:)/N,'s-','Color',ccc(1,:),'LineWidth',2 )
plot(1:1:Hours,average_consumption_24(4,:)/N,'s-','Color',ccc(3,:),'LineWidth',2 )
plot(1:1:Hours,average_consumption_24(6,:)/N,'s-','Color',ccc(5,:),'LineWidth',2 )
plot(1:1:Hours,average_consumption_24(8,:)/N,'s-','Color',ccc(7,:),'LineWidth',2 )
grid on
title('Hybrid Models - 24 Hours')
ylabel('Average Consumption per Capita')
xlabel('Hours')
lll = legend('$\Gamma = $H7, $\Omega = $B','$\Gamma = $H8, $\Omega = $B',...
    '$\Gamma = $H9, $\Omega = $B','$\Gamma = $H10, $\Omega = $B');
set(lll,'Interpreter','Latex');
set(lll,'FontSize',14);
hold off
