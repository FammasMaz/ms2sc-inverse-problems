n1 = 0.5;
n2 = 1.5;

elast = dsigma*t_range/E_ref;
platic1 = ((dsigma/K_ref)^n1)*(t_range.^(n1+1))/(n1+1);
platic2 = ((dsigma/K_ref)^n2)*(t_range.^(n2+1))/(n2+1);

figure(1)
plot(t_range,elast,t_range,platic1,t_range,platic2)
legend('elast','n < 1','n > 1')