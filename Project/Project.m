clear; close all; clc;

% Dataset
% Define parameters (by default, these are the initial guesses for optimizations)
beam.E = 300e+09; % Pa
beam.K = 90e+06; % Pa
beam.n = 12;
dsigma = 0.8e6; % Pa*s^-1
r = 0.0254; % m
s = pi * r^2; % m^2

% Set range of parameters
n_range = [4:0.1:30];
e_range = [250:10:400].*1e+09;
k_range = [70:1:120].*1e+06;

% Define time span
tspan = [0:0.1:80]'; % s

% Noise level in percentage
noise_level = 0;

% Define ODE function
%[t, epsilon_exp]=ode45(@forward_sigma, tspan, 0, [], beam, dsigma);
epsilon_exp = forana(beam, tspan, dsigma); % creation of synthetic data
t = tspan;
sigma_exp = dsigma * t; % creation of synthetic data
N = length(epsilon_exp);


% Plot the forward problem
figure
plot(epsilon_exp, sigma_exp);
if noise_level ~= 0
    hold on;
    epsilon_exp = epsilon_exp + noise_level/100 * randn(N,1);
    plot(epsilon_exp, sigma_exp);
    hold off;
    legend('Original', 'Noisy');
end
xlabel('Strain');
ylabel('Stress');
title('Stress-Strain Curve: Forward Problem with Guessed Parameters');
saveas(gcf, 'assets/forward_problem.png');


% Plot misfit values for different values of n
misfit_plotter(e_range, k_range, n_range, epsilon_exp, t, beam, dsigma, 'n');

% Plot misfit values for different values of E
misfit_plotter(e_range, k_range, n_range, epsilon_exp, t, beam, dsigma, 'E');

% Plot misfit values for different values of K
misfit_plotter(e_range, k_range, n_range, epsilon_exp, t, beam, dsigma, 'K');


% Global optimization using manually iterating over the range of parameters
[misfit_values, E_opt, K_opt, n_opt] = misfit_global(epsilon_exp, t, beam, dsigma, e_range, k_range, n_range);
% Define objective function
fun = @(x) misfit_sig(x, epsilon_exp, t, beam, dsigma);

% Define initial guess and search bounds
x0 = [100e9; 60e7; 12]; % initial guess

% Optimization
options = optimset('Display','iter'); % display iterations
[x_opt, fval] = fminsearch(fun, x0, options);   % local optimization

% Display results
fprintf('Optimal parameter values using fminsearch:\n');
fprintf('E = %g Pa\n', x_opt(1));
fprintf('K = %g Pa\n', x_opt(2));
fprintf('n = %g\n', x_opt(3));
fprintf('Misfit function value = %g\n', fval);

% Generate mesh of E and n ranges
mesh_plotter(misfit_values, e_range, k_range, n_range, 'K', beam);

% Generate mesh of K and n ranges
mesh_plotter(misfit_values, e_range, k_range, n_range, 'E', beam);

% Generate mesh of K and E ranges
mesh_plotter(misfit_values, e_range, k_range, n_range, 'n', beam);


%% Inverse Problem (Identification usign fminsearch)

% fminsearch used since from the shape of misfit suggests that a global
% gradient based optimization is not suitable since the optimizer can get
% stuck in a local minima. `fminsearch` follows a downhill simplex method 
% introduced by Nelder and Mead in 1965 (IP-Session3, slide 26 on Course notes)

% Load id1.mat (two variables: epsiexp and sigexp)
load('id1.mat');

% Define objective function
fun = @(x) misfit_sig(x, epsiexp, t, beam, dsigma);

% Define initial guess and search bounds
x0 = [beam.E; beam.K; beam.n]; % initial guess

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define objective function
%%
% Run optimization


%% Other tests that didnt work

% It is found that the misfit function is not convex. Hence, the global optimization is not possible.
% The purpose of using fminsearch instead of the fminunc is to find a solution without checking for gradients

% %[misfit, E_opt, K_opt, n_opt] = misfit_minunc(epsilon_exp, t, beam, dsigma, e_range, k_range, n_range);
% % e_range = 180:10:300;
% % e_range = e_range .* 1e9;
% % k_range = 65:1:75;
% % k_range = k_range .* 1e6;
% %[misfit_values, E_opt, K_opt, n_opt] = misfit_mincon(epsilon_exp, t, beam, dsigma, e_range, [beam.K beam.K], [beam.n, beam.n]);


% %id1 mat
% tid = [0:0.1:80];

% load('id1.mat');
% E_guess = 200e9;
% K_guess = 70e7;
% n_guess = 3;
% x0 = [E_guess, K_guess, n_guess];

% %[misfit, E_opt, K_opt, n_opt] = ident(x0, epsiexp, tid, beam, dsigma);

% tn = [0:0.1:80]';
% epsilonid = forana(beam, tn, dsigma);
% sigmaid = dsigma * tn;

% [xsol,fval]=fminunc(@misfitid, x0,optimset('Display','iter','TolFun',1e-6,'GradObj','off'),epsiexp, tn, beam, dsigma) % unconstrained minimization of the misfit function
% [misfit_values, E_opt, K_opt, n_opt] = misfit_mincon(epsiexp, tn, beam, dsigma, [E_guess], [K_guess], [n_guess]);


% % Define experimental data
% load('id1.mat');

% tn = 0:0.1:80;
% dsigma = 0.8e6 % your value of dsigma here

% % Set initial guesses for parameters
% initial_guess = [200e9, 70e6, 1.5];

% % Set up the options for the optimization algorithm
% options = optimoptions('fminunc', 'MaxIterations', 1000);

% % Minimize the misfit function using fminunc
% result = fminunc(@(beam) misfitid(beam, tn, epsiexp, dsigma), initial_guess, options);