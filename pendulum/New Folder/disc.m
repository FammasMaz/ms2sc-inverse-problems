function disc = disc(t, theta, thetaexp, texp)
    ind = dsearchn(t, texp');

    disc = 0.5*(theta(ind,1) - thetaexp(ind,1)).^2;

