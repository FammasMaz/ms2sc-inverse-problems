

% Define the time interval
t = linspace(0, 10, 1000);

% Call the function to solve the ODE
solution = solve_de(t);

% Plot the solution
plot(solution(:,1), solution(:,2));
xlabel('Time (s)');
ylabel('Strain rate (s^{-1})');


function solution = solve_de(t)
    
    E = 200000; % Modulus of elasticity
    K = 1000; % Material constant
    n = 0.2; % Material constant
    
    r = 2.54; % radius in cm
    s = pi * r^2; % cross-sectional area
    
    dot_sigma = 0.8 / s; % stress rate
    
    % Define the right-hand side function of the ODE
    % and the initial condition
    
    ode_func = @(t, epsilon) (dot_sigma/E) + ((sigma(t)/K)^n);
    ic = [0];
    
    [t, epsilon] = ode45(ode_func, t, ic);
    
    solution = [t, epsilon];
end