% Parameter identification: E, K, n

% Information
dsigma = 0.8e6;

% Experimental data
sigma_data = sigexp';
epsilon_data = epsiexp';

t_final = sigma_data(end)/dsigma;
t_step = t_final/(length(sigma_data)-1);
t_range = 0:t_step:t_final;

% Minimization of the misfit function
x_init = [200e9 100e6 8];
sol = fminunc(@misfit,x_init,optimset('Display','iter','TolFun',1e-6,'GradObj','off'),dsigma,t_range,epsilon_data);
E_predict = sol(1);
K_predict = sol(2);
n_predict = sol(3);

sigma_predict = fwd_sigma(t_range,dsigma);
x_data = [E_predict K_predict n_predict];
epsilon_predict = fwd_epsilon(x_data,t_range,dsigma);

figure(1);
plot(epsilon_data,sigma_data,epsilon_predict,sigma_predict);
xlabel("strain");
ylabel("stress");
legend('experiment','prediction')
