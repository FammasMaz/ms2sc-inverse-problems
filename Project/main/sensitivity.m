function [sens, param] = sensitivity(x0, noise_level, epsilon_exp, t, beam, dsigma, lambda, arg)

options = optimset(); % display iterations
sens = zeros(length(noise_level), 1);
param = zeros(length(noise_level), length(x0));
iter = zeros(length(noise_level), 1);
for i = 1:length(noise_level)
    epsilon_exp_noise = epsilon_exp + noise_level(i)*randn(size(epsilon_exp));
    fun = @(x) misfit_sig(x, epsilon_exp_noise, t, beam, dsigma, lambda, arg);
    [x_opt, fval, exitflag, output] = fminsearch(fun, x0, options); 
    sens(i) = fval;
    param(i,:) = x_opt;
    iter(i) = output.iterations;
end

figure
plot(noise_level, param(:,1), 'r');
hold on;
plot(noise_level, param(1,1)*ones(size(noise_level)), 'b');
hold off;
xlabel('Noise level (%)');
ylabel('E');
legend('Perturbed', 'Original');
title('fminsearch: Sensitivity of E to noise');
saveas(gcf, ['assets/sensitivity_E',arg,'.png']);

figure
plot(noise_level, param(:,2), 'r');
hold on;
plot(noise_level, param(1,2)*ones(size(noise_level)), 'b');
hold off;
xlabel('Noise level (%)');
ylabel('K');
legend('Perturbed', 'Original');
title('fminsearch: Sensitivity of K to noise');
saveas(gcf, ['assets/sensitivity_K',arg,'.png']);

figure
plot(noise_level, param(:,3), 'r');
hold on;
plot(noise_level, param(1,3)*ones(size(noise_level)), 'b');
hold off;
xlabel('Noise level (%)');
ylabel('n');
legend('Perturbed', 'Original');
title('fminsearch: Sensitivity of n to noise');
saveas(gcf, ['assets/sensitivity_n',arg,'.png']);

figure
plot(noise_level, iter, 'r');
xlabel('Noise level (%)');
ylabel('Iterations');
title('fminsearch: Iterations to converge');
saveas(gcf, ['assets/sensitivity_iter',arg,'.png']);

end


