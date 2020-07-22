clear
clc
close all

%% Initilization
N = 100;
runtime = 50;

%% Graph
count = 0;
maxnum = 0;
minnum = 100;
for run = 1:runtime
    graph = generate_graph(N,run);
    incidence_matrix = graph;
    nr_nodes = N;
    nr_neighbors = incidence_matrix*ones(nr_nodes,1);
    count = count + mean(nr_neighbors);
    if max(nr_neighbors) > maxnum
        maxnum = max(nr_neighbors);
    end
    if min(nr_neighbors) < minnum
        minnum = min(nr_neighbors);
    end
end

disp(count/runtime)
    