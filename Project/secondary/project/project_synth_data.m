% Parameter identification: E, K, n
clear all; close all; clc;

% Information
dsigma = 0.8e6;

% Parameters
E_range = 250e9:10e9:400e9;
K_range = 70e6:1e6:120e6;
n_range = 4:0.1:30;

E_data = 300e9;
K_data = 90e6;
n_data = 12;

t_range = 0:1:100;

% Synthetic data
sigma_data = fwd_sigma(t_range,dsigma);
x_data = [E_data K_data n_data];
epsilon_data = fwd_epsilon(x_data,t_range,dsigma);

% Noisy synthetic data
% epsilon_data = epsilon_data + 10^(-1)*randn(1,length(t_range));

% Parameter sweep to plot the misfit function
% matrix_misfit = sweep_misfit(E_range(1),K_range,n_range,dsigma,t_range,epsilon_data);
% matrix_misfit = reshape(matrix_misfit,length(K_range),length(n_range));
% 
% figure(2)
% [X_Kn,Y_Kn] = meshgrid(K_range,n_range);
% surf(X_Kn,Y_Kn,matrix_misfit');
% xlabel('K'); ylabel('n'); zlabel('misfit');

% Minimization of the misfit function
% x_init = [250e9 70e6 9];

% sol_unc = fminunc(@misfit,x_init,optimset('Display','iter','TolFun',1e-6,'GradObj','off'),dsigma,t_range,epsilon_data);
% E_unc = sol_unc(1);
% K_unc = sol_unc(2);
% n_unc = sol_unc(3);
% sigma_unc = fwd_sigma(t_range,dsigma);
% x_unc = [E_unc K_unc n_unc];
% epsilon_unc = fwd_epsilon(x_unc,t_range,dsigma);

% figure(1);
% plot(epsilon_data,sigma_data,epsilon_unc,sigma_unc,'+');
% xlabel("strain");
% ylabel("stress");
% legend('synthetic experiment','prediction fminunc')
