close all; clear all; clc;
% Parameter identification

% Experimental data, parameters
E = 200e9; % Pa
K = 5e6; % Pa
n = 2.75;
sigma_dot = 0.8e6; % Pa*s^-1
r = 0.0254; % m
s = pi * r^2; % m^2

% ode45 uses runge kutta order 4
[t, strainexp]=ode45(@forward,[0:.01:1],0,[],K, n, s, E, sigma_dot); % creation of synthetic data
figure;plot(t,strainexp,'r')
load("id1.mat");

%strainexp = epsiexp;
% Initial model
n0=1; % initial parameter value for the forward model
[t, strain]=ode45(@forward,[0:.01:1],0,[],K, n0, s, E, sigma_dot); % creation of synthetic data
hold on;plot(t,strain, 'g')
stress = sigma_dot*t
figure;plot(strainexp, stress,"g")
% Parameter sweep to plot the misfit function
for i=1:30
    value(i)=misfit(i/10,strainexp,t, K, s, E, sigma_dot);
end
figure;plot([.1:.1:3],value)

% Minimization of the misfit function
[nsol,fval]=fminunc(@misfit,n0,optimset('Display','iter','TolFun',1e-6,'GradObj','off'),strainexp,t, K, s, E, sigma_dot) % unconstrained minimization of the misfit function
[tsol,strainsol]=ode45(@forward,[0:.01:10], 0,[], K, nsol, s, E, sigma_dot); % predictions of the identified model
figure;plot(t,strainexp,'r')
hold on;plot(tsol,strainsol)
