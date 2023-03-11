clear all; close all; clc;

%% Forward Problem
% Define parameters
beam.E = 200e9; % Pa
beam.K = 75e6; % Pa
beam.n = 30.8560;
dsigma = 0.8e6; % Pa*s^-1
r = 0.0254; % m
s = pi * r^2; % m^2

% Define time span
tspan = [0:0.1:80]; % s

% Define ODE function
[t, epsilon_exp]=ode45(@forward_sigma, tspan, 0, [], beam, dsigma); % creation of synthetic data
sigma_exp = dsigma * t; % creation of synthetic data
N = length(epsilon_exp);
%epsilon_exp = epsilon_exp + 0.01 * randn(N,1);


% Plot the forward problem
figure
plot(epsilon_exp, sigma_exp);
xlabel('Strain');
ylabel('Stress');
title('Stress-Strain Curve: Forward Problem');

n_range = 0.5:0.1:2.5;
for i=1:length(n_range)   
    x = [beam.E; beam.K; n_range(i)]; 
    value(i)=misfit_sig(x, epsilon_exp, t, beam, dsigma);
end

%% Misfit Function Under Testing
% x0 = [270e9; 72e6; 0.13];
% obj_func = @(x) misfit_sig(x, epsilon_exp, t, beam, dsigma);
% options = optimset('Display','iter','TolFun',1e-6,'GradObj','off');
% [x_opt, misfit_min] = fminunc(obj_func, x0, options);
%%%




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
mesh(E_mesh, n_mesh, squeeze(misfit_values(:,end,:))); % plot for n = n_range(1)
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


%[misfit, E_opt, K_opt, n_opt] = misfit_minunc(epsilon_exp, t, beam, dsigma, e_range, k_range, n_range);
% e_range = 180:10:300;
% e_range = e_range .* 1e9;
% k_range = 65:1:75;
% k_range = k_range .* 1e6;
%[misfit_values, E_opt, K_opt, n_opt] = misfit_mincon(epsilon_exp, t, beam, dsigma, e_range, [beam.K beam.K], [beam.n, beam.n]);


%id1 mat
tid = [0:0.1:80];

load('id1.mat');
E_guess = 200e9;
K_guess = 70e7;
n_guess = 3;
x0 = [E_guess, K_guess, n_guess];

%[misfit, E_opt, K_opt, n_opt] = ident(x0, epsiexp, tid, beam, dsigma);

tn = [0:0.1:80]';
epsilonid = forana(beam, tn, dsigma);
sigmaid = dsigma * tn;

[xsol,fval]=fminunc(@misfitid, x0,optimset('Display','iter','TolFun',1e-6,'GradObj','off'),epsiexp, tn, beam, dsigma) % unconstrained minimization of the misfit function
[misfit_values, E_opt, K_opt, n_opt] = misfit_mincon(epsiexp, tn, beam, dsigma, [E_guess], [K_guess], [n_guess]);
