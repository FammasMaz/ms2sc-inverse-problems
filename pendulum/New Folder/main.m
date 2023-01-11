
% Experimental Data
gexp = 10;
[texp, thetaexp] = ode45(@forw, [0:0.1:10], [pi/3 0], [], gexp);

plot(texp, thetaexp(:,1),'r');

hold on;

% Initial Data
g = 7
[tinit, thetainit] = ode45(@forw, [0:0.1:10], [pi/3 0], [], g);

plot(tinit, thetainit(:,1),'b');

hold off;


for i=1:200
    value(i)=misf(i/10,thetaexp,texp);
end
figure;plot([.1:.1:20],value)
