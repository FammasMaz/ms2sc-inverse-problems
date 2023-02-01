function forward=forward(t,F, K, n, s, E, sigma_dot)
% % l=1 : a unit length is assumed here
% forward=[theta(2);-g*sin(theta(1))]; % in the form y'=f(y,t)

forward = sigma_dot / E + (sigma_dot *t / K)^n;

