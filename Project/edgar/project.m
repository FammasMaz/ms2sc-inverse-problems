% Parameter identification: E, K, n
clear all; close all; clc;

% Information
dsigma = 0.8e6;

% Parameters
E_range = 100e9:10e9:300e9;
K_range = 10e6:5e6:90e6;
n_range = 0.1:0.1:2;

E_ref = 0.5*(E_range(1)+E_range(end));
K_ref = 0.5*(K_range(1)+K_range(end));
n_ref = 0.5*(n_range(1)+n_range(end));

E_data = 180e9;
K_data = 70e6;
n_data = 2;
t_range = 0:0.01:10;

% Synthetic data
sigma_data = fwd_sigma(t_range,dsigma);
epsilon_data = fwd_epsilon(t_range,dsigma,E_data,K_data,n_data);

% figure(1);
% plot(epsilon_data,sigma_data);
% xlabel("strain");
% ylabel("stress");

% Parameter sweep to plot the misfit function
val_misfit = misfit(dsigma,t_range,epsilon_data,E_range,K_range,n_ref);

[X,Y] = meshgrid(E_range,K_range);
surf(X,Y,val_misfit');
