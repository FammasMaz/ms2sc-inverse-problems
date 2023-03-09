function forward=forward(t, beam)
% l=1 : a unit length is assumed here
forward=[theta(2);-g*sin(theta(1))]; % in the form y'=f(y,t)


forward = [((dsigma * t / beam.K)^beam.n) + (0.8/beam.E)]