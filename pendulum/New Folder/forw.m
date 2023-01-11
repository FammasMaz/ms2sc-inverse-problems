function forw = forw(t, theta, g)
    forw = [theta(2); -g*sin(theta(1))/1];
end
