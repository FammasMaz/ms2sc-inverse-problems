function misfitid=misfitid(x, eps_exp, tn, beam, dsigma)
% Calculation of the misfit function with no gradient estimation
beamtemp = beam;
beamtemp.n = x(3); % update of the parameter
beamtemp.E = x(1);
beamtemp.K = x(2);
beamtemp
epsilonid = forana(beamtemp, tn, dsigma);
misfitid=quad(@discrepancy,tn(1),tn(end),[],0,eps_exp,epsilonid,x,tn);