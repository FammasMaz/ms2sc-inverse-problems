function misfit=misfit(g,thetaexp,tn)
% Calculation of the misfit function with no gradient estimation
[t,theta]=ode45(@forward,tn,[pi/3 0],[],g); % simulation of the response using the parameter g
misfit=quad(@discrepancy,tn(1),tn(end),[],0,thetaexp,theta,g,tn); % calculation of the integral over time of the discrepancy model-exp