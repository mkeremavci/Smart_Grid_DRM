function [ graph ] = generate_graph(N,run)
    graph = zeros(N, N);
    rng(run)
    x_axis = 8*rand(1, 100);
    rng(run)
    y_axis = 8*rand(1, 100);
    point = x_axis +1j.*y_axis;
    for i=1:N
       for j=1:N
           if abs(point(i)-point(j))<2
              graph(i,j) = 1; 
           end
       end
    end
end