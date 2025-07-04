function f = initialize_variables(N, M, V, min_range, max_range)
lb = min_range;
ub = max_range;
K = M + V;

for i = 1 : N
    for j = 1 : V
        f(i,j) = lb(j) + (ub(j) - lb(j))*rand(1);
                                                   
    end
    f(i,V + 1: K) = evaluate_objective(f(i,:), M, V);
                                                   
end