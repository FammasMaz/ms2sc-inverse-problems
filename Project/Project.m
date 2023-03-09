clear all; close all; clc;
% Define parameters
beam.E = 200e9; % Pa
beam.K = 70e6; % Pa
beam.n = 2;
dsigma = 0.8e6; % Pa*s^-1
r = 0.0254; % m
s = pi * r^2; % m^2

% Define time span and step size
tspan = [0, 10]; % s
step = 0.01; % s

% Define ODE function
% [t, eps] = ode45(@(t, eps) ((dsigma * t / beam.K)^beam.n)+(0.8/beam.E),[0:.01:10], 0);
[t, epsilon]=ode45(@forward_sigma, [0:.01:10], 0, [], beam, dsigma); % creation of synthetic data

figure
plot(epsilon, sigma);
xlabel('Strain');
ylabel('Stress');
title('Stress-Strain Curve: Forward Problem');