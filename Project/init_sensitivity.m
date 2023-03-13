function [iter, wrong] = init_sensitivity(e_range, k_range, n_range, fun, algo, beam, arg)
    iter = zeros(length(e_range), length(k_range), length(n_range));
    wrong = 0;
    for i = 1:length(e_range)
        for j = 1:length(k_range)
            for k = 1:length(n_range)
                x0 = [e_range(i), k_range(j), n_range(k)];
                [x, ~, ~, output] = fminsearch(fun, x0, optimset());
                iter(i, j, k) = output.iterations;
                test = isequal(round(x(3)), round(beam.n));
                if test == 0
                    disp('Solution guess not equal to initial guess');
                    iter(i, j, k) = 500;
                    wrong = wrong+1;
                end
            end
        end
    end
    figure
    surf(e_range, k_range, squeeze(iter(:, :, 1))')
    xlabel('E')
    ylabel('K')
    zlabel('Iterations')
    title('Iterations for different E and K values')
    saveas(gcf, strcat('assets/init_sens', '_', algo, '_e_k',arg,'.png'))
    figure
    surf(e_range, n_range, squeeze(iter(:, 1, :))')
    xlabel('E')
    ylabel('n')
    zlabel('Iterations')
    title('Iterations for different E and n values')
    saveas(gcf, strcat('assets/init_sens', '_', algo, '_e_n',arg,'.png'))
    figure
    surf(k_range, n_range, squeeze(iter(1, :, :))')
    xlabel('K')
    ylabel('n')
    zlabel('Iterations')
    title('Iterations for different K and n values')
    saveas(gcf, strcat('assets/init_sens', '_', algo, '_k_n', arg,'.png'))

    fprintf('Wrong solutions: %d', wrong);
end

