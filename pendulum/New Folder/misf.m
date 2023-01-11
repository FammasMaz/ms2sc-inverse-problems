function misf=misf(g, thetaexp, texp)

    [t, theta] = ode45(@forw, texp, [pi/3 0], [], g);

    misf = quad(@disc,texp(1),texp(end),[],0,theta,thetaexp,texp');
end