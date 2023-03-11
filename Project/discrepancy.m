function discrepancy=discrepancy(t,thetaexp,theta,g,tn)
% calculation of the misfit function at time t
ind=dsearchn(tn,t'); % ind is the vector of the indices of time steps in tn that are the closest to the time steps contained in the vector t
% alternative to the function dsearchn below
%for i=1:length(t)
%[m,ind(i)]=min(abs(t(i)-tn));
%end
discrepancy=0.5*(theta(ind,1)-thetaexp(ind,1)).^2; % quantity to be integrated over time
