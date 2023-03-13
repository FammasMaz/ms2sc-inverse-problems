clear; close all; clc;

% Dataset
% Define parameters (by default, these are the initial guesses for optimizations)
beam.E = 300e+09; % Pa
beam.K = 90e+06; % Pa.s
beam.n = 12;
dsigma = 0.8e6; % Pa*s^-1
r = 0.0254; % m
s = pi * r^2; % m^2
arg = 'RMSE'

% Set range of parameters
n_range = [4:0.5:30];
e_range = [250:10:400].*1e+09;
k_range = [70:1:120].*1e+06;
xtest = [50e9, 100e6, 9]; % initial guess
x0 = [100e9; 60e6; 11]; % initial guess

% Define time span
tspan = [0:0.1:80]'; % s
t = tspan;
% Noise level in percentage and regularization
noise_level = 0;
lambda = 0;
% Define ODE function
%[t, epsilon_exp]=ode45(@forward_sigma, tspan, 0, [], beam, dsigma);
epsilon_exp = forana(beam, tspan, dsigma); % creation of synthetic data
t = tspan;
sigma_exp = dsigma * t; % creation of synthetic data
N = length(epsilon_exp);


% Plot the forward problem
figure
plot(epsilon_exp, sigma_exp);
xlim([-0.02, 0.15]);
tit = 'assets/forward_problem.png';
if noise_level ~= 0
    hold on;
    epsilon_exp = epsilon_exp + noise_level/100 * randn(N,1);
    plot(epsilon_exp, sigma_exp);
    if lambda~=0
        eps_fil = low_pass_filter(tspan, epsilon_exp, lambda);
        plot(eps_fil, sigma_exp, 'g');
        legend('Original', 'Noisy', 'Filtered');
    else
        legend('Original', 'Noisy');
    end
    hold off;
    tit = 'assets/forward_problem_noisy.png';
end
xlabel('Strain');
ylabel('Stress');
title('Stress-Strain Curve: Forward Problem with Guessed Parameters');
saveas(gcf, tit);

%% Effect of the time step
time_effect(x0, beam, dsigma, lambda, arg)


%% Sensitivity Analysis against noise
noise_range = [0.0:0.05:2];
%[sens, param] = sensitivity(xtest, noise_range, epsilon_exp, t, beam, dsigma, lambda, arg);


%% Misfit Plotter
% Plot misfit values for different values of n
misfit_plotter(e_range, k_range, n_range, epsilon_exp, t, beam, dsigma, 'n', lambda,arg);
% Plot misfit values for different values of E
misfit_plotter(e_range, k_range, n_range, epsilon_exp, t, beam, dsigma, 'E', lambda,arg);
% Plot misfit values for different values of K
misfit_plotter(e_range, k_range, n_range, epsilon_exp, t, beam, dsigma, 'K',lambda, arg);

%% Global optimization using misfit_global
% Global optimization using manually iterating over the range of parameters
%[misfit_values, E_opt, K_opt, n_opt] = misfit_global(epsilon_exp, t, beam, dsigma, e_range, k_range, n_range, lambda, arg);

%% Sensitivity Analysis to initial guess
% Define initial guess and search bounds

% Plotting the iterations taken depending on the initial guess
options = optimset('Display','iter'); % display iterations
fun = @(x) misfit_sig(x, epsilon_exp, t, beam, dsigma, lambda, arg);
e_init = [200e9:10e9:230e9];
k_init = [60e6:10e6:80e6];
n_init = [10:1:15];
%[iter, wrong] = init_sensitivity(e_init, k_init, n_init, fun, 'fminsearch', beam, arg);

% Optimization
fun = @(x) misfit_sig(x, epsilon_exp, t, beam, dsigma, lambda, arg);
[x_opt, fval] = fminsearch(fun, x0, options); 
[xsol, fval]=fminunc(fun, x0, optimset('Display','iter','TolFun',1e-6,'GradObj','off'));

% Display results
fprintf('Optimal parameter values using fminsearch:\n');
fprintf('E = %g Pa\n', x_opt(1));
fprintf('K = %g Pa\n', x_opt(2));
fprintf('n = %g\n', x_opt(3));
fprintf('Misfit function value = %g\n', fval);

%% Misfit plotter 3D meshes
% Generate mesh of E and n ranges
mesh_plotter(misfit_values, e_range, k_range, n_range, 'K', beam, arg);
% Generate mesh of K and n ranges
mesh_plotter(misfit_values, e_range, k_range, n_range, 'E', beam, arg);
% Generate mesh of K and E ranges
mesh_plotter(misfit_values, e_range, k_range, n_range, 'n', beam, arg);


%% Inverse Problem (Identification usign fminsearch)

% fminsearch used since from the shape of misfit suggests that a global
% gradient based optimization is not suitable since the optimizer can get
% stuck in a local minima. `fminsearch` follows a downhill simplex method 
% introduced by Nelder and Mead in 1965 (IP-Session3, slide 26 on Course notes)

% Load id1.mat (two variables: epsiexp and sigexp)
load('id1.mat');

% Define objective function
fun = @(x) misfit_sig(x, epsiexp, t, beam, dsigma, arg);

% Define initial guess and search bounds
x0 = [beam.E; beam.K; beam.n]; % initial guess
x0 = [300e9; 80e6; 12];

% Optimization
options = optimset('Display','iter'); % display iterations
[x_opt, fval] = fminsearch(fun, x0, options);   % local optimization

% Display results
fprintf('Optimal parameter values:\n');
fprintf('E = %g Pa\n', x_opt(1));
fprintf('K = %g Pa\n', x_opt(2));
fprintf('n = %g\n', x_opt(3));
fprintf('Misfit function value = %g\n', fval);

%% Plot using optimal parameters
idbeam.E = x_opt(1);
idbeam.K = x_opt(2);
idbeam.n = x_opt(3);
epsilon_id = forana(idbeam, t, dsigma); % creation of synthetic data using the identified model
sigma_id = dsigma * t; % creation of synthetic data using the identified model

figure
plot(epsiexp, sigexp, 'b', epsilon_id, sigma_id, 'r');
xlabel('Strain');
ylabel('Stress');
title('Stress-Strain Curve: Data vs Identified Model');
legend('Data', 'Identified Model');
saveas(gcf, 'assets/identified_model.png');