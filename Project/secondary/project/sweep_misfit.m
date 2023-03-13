function matrix_misfit = sweep_misfit(E,K,n,dsigma,t_range,epsilon_data)

    matrix_misfit = zeros(length(E),length(K),length(n));
    for i = 1:length(E)
        for j = 1:length(K)
            for h = 1:length(n)
                x_range = [E(i) K(j) n(h)];
                matrix_misfit(i,j,h) = misfit(x_range,dsigma,t_range,epsilon_data);
            end
        end
    end

end