function discr = discrepancy (t,x,dsigma,t_range,epsilon_data)
    eps = fwd_epsilon(x,t,dsigma);
    ind=dsearchn(t_range',t');
    discr = 0.5*(eps - epsilon_data(ind)).^2;
end