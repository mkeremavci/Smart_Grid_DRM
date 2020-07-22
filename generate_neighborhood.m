function [ graph ] = generate_neighborhood(N,num)
    graph = zeros(N, N);
    axis = zeros(1,num);
    axis(1:mod(N,num)) = ceil(N/num);
    axis = axis + (axis == 0) .* floor(N/num);
    axis = cumsum(axis);
    for i = 1:num
        if i == 1
            graph(1:axis(i),1:axis(i)) = 1;
            continue
        end
        graph(axis(i-1)+1:axis(i),axis(i-1)+1:axis(i)) = 1;
    end
end