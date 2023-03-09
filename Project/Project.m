clear all; close all; clc;

%% Forward Problem
% Define parameters
beam.E = 200e9; % Pa
beam.K = 70e6; % Pa
beam.n = 1.1;
dsigma = 0.8e6; % Pa*s^-1
r = 0.0254; % m
s = pi * r^2; % m^2

% Define time span
tspan = [0:0.01:10]; % s

% Define ODE function
[t, epsilon_exp]=ode45(@forward_sigma, tspan, 0, [], beam, dsigma); % creation of synthetic data
sigma_exp = dsigma * t; % creation of synthetic data

% Plot the forward problem
figure
plot(epsilon_exp, sigma_exp);
xlabel('Strain');
ylabel('Stress');
title('Stress-Strain Curve: Forward Problem');

n_range = 0.5:0.1:2.5;

for i=1:length(n_range)    
    value(i)=misfit_sig(n_range(i), epsilon_exp, t, beam, dsigma);
end

figure
plot(n_range, value);
e_range = 50:10:250;
e_range = e_range .* 1e9;
k_range = 65:1:75;
k_range = k_range .* 1e6;

[misfit_values, E_opt, K_opt, n_opt] = misfit_global(epsilon_exp, t, beam, dsigma, e_range, k_range, n_range);

% Generate mesh of E and n ranges
[E_mesh, n_mesh] = meshgrid(e_range, n_range);

% Create 3D mesh plot of misfit values
figure
mesh(E_mesh, n_mesh, squeeze(misfit_values(:,5,:))); % plot for n = n_range(1)
xlabel('E');
ylabel('n');
zlabel('Misfit');
title('Misfit Values for Different Values of E and n');

% Generate mesh of K and n ranges
[K_mesh, n_mesh] = meshgrid(k_range, n_range);

% Create 3D mesh plot of misfit values
figure
mesh(K_mesh, n_mesh, squeeze(misfit_values(1,:,:))'); % plot for n = n_range(1)
xlabel('K');
ylabel('n');
zlabel('Misfit');
title('Misfit Values for Different Values of K and n');


% Generate mesh of K and E ranges
[K_mesh, E_mesh] = meshgrid(k_range, e_range);
% Create 3D mesh plot of misfit values
figure
mesh(K_mesh, E_mesh, squeeze(misfit_values(:,:,5))); % plot for n = n_range(1)
xlabel('K');
ylabel('E');
zlabel('Misfit');
title('Misfit Values for Different Values of K and E');
