clear all; close all; clc;


addpath("misfit/");
load('id1.mat');

plot(epsiexp, sigexp);
xlim([0, 0.1])


function forward=forward(t,theta,g)
% l=1 : a unit length is assumed here
forward=[theta(2);-g*sin(theta(1))]; % in the form y'=f(y,t)


% Define parameters
E = 100e9; % Pa
K = 1e6; % Pa
n = 0.2;
sigma_dot = 0.8e6; % Pa*s^-1
r = 0.0254; % m
s = pi * r^2; % m^2

% Define time span and step size
tspan = [0, 10]; % s
step = 0.01; % s

% Define ODE function
ode_fun = @(t, F) sigma_dot / E + (F / (s * K))^n;

% Solve ODE
[t, F] = ode45(ode_fun, tspan, 0);

% Plot the result
plot(t, F)
xlabel('time (s)')
ylabel('F(t) (N)')

