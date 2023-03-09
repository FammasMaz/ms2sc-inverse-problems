function forward_sigma=forward_sigma(t, eps, beam, dsigma)
% l=1 : a unit length is assumed here
forward_sigma = [((dsigma * t / beam.K)^beam.n)+(0.8/beam.E)];