function forward=forward(t,theta,g)
% l=1 : a unit length is assumed here
forward=[theta(2);-g*sin(theta(1))]; % in the form y'=f(y,t)


