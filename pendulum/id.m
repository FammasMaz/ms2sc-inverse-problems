% Parameter identification

% Experimental data
gexp=10; % `true’ value of the parameter g

% ode45 uses runge kutta order 4
[t,thetaexp]=ode45(@forward,[0:.01:10],[pi/3 0],[],gexp); % creation of synthetic data
figure;plot(t,thetaexp(:,1),'r')

% Initial model
g0=8; % initial parameter value for the forward model
[t,theta]=ode45(@forward,[0:.01:10],[pi/3 0],[],g0); % predictions of the model
hold on;plot(t,theta(:,1))

% Parameter sweep to plot the misfit function
for i=1:200
    value(i)=misfit(i/10,thetaexp,t);
end
figure;plot([.1:.1:20],value)

% Minimization of the misfit function
[gsol,fval]=fminunc(@misfit,g0,optimset('Display','iter','TolFun',1e-6,'GradObj','off'),thetaexp,t) % unconstrained minimization of the misfit function
[t,thetasol]=ode45(@forward,[0:.01:10],[pi/3 0],[],gsol); % predictions of the identified model
figure;plot(t,thetaexp(:,1),'r')
hold on;plot(t,thetasol(:,1))
