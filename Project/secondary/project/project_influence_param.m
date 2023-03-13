% Parameter identification: E, K, n
clear all; close all; clc;

% Information
dsigma = 0.8e6;

% Parameters
E_init = 50e9;
K_init = 50e6;
n_init = 2;

n_step = 0.1;
E_step = n_step*1e9;
K_step = n_step*1e6;

t = 50;

% Synthetic data
sigma_data = fwd_sigma(t,dsigma);
x_data = [E_init K_init n_init];
epsilon_data_init = fwd_epsilon(x_data,t,dsigma);

% Morris' method
Q = 200;
mean_effect = zeros(1,3);
deviation_effect = zeros(1,3);

elem_effect = zeros(1,Q);
for i = 1:Q
    x_data = [E_init+(i*E_step) K_init n_init];
    epsilon_data = fwd_epsilon(x_data,t,dsigma);
    elem_effect(i) = (epsilon_data - epsilon_data_init)/(i*n_step);
end
mean_effect(1) = mean(elem_effect);
deviation_effect(1) = sqrt(mean((elem_effect-mean_effect(1)).^2));

elem_effect = zeros(1,Q);
for i = 1:Q
    x_data = [E_init K_init+(i*K_step) n_init];
    epsilon_data = fwd_epsilon(x_data,t,dsigma);
    elem_effect(i) = (epsilon_data - epsilon_data_init)/(i*n_step);
end
mean_effect(2) = mean(elem_effect);
deviation_effect(2) = sqrt(mean((elem_effect-mean_effect(2)).^2));

elem_effect = zeros(1,Q);
for i = 1:Q
    x_data = [E_init K_init n_init+(i*n_step)];
    epsilon_data = fwd_epsilon(x_data,t,dsigma);
    elem_effect(i) = (epsilon_data - epsilon_data_init)/(i*n_step);
end
mean_effect(3) = mean(elem_effect);
deviation_effect(3) = sqrt(mean((elem_effect-mean_effect(3)).^2));

plot(abs(mean_effect(1)),deviation_effect(1),'+',abs(mean_effect(2)),deviation_effect(2),'+',abs(mean_effect(3)),deviation_effect(3),'+')
xlabel('mean')
ylabel('standard deviation')
legend('E','K','n')