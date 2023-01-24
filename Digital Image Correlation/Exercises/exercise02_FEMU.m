clear; close all

% load some functions provided in the lib folder
addpath('lib');

% Exercise 2 : FEMU
% -----------------------------------------------------------

% FEMU is short for Finite Element Method Updating. It is an inverse method
% where we iteratively optimise the material parameters until we have
% minimise the distance with the experiment.

% The result of the previous exercise is the displacement field on the
% surface of the sample. Here it is loaded from file, together with the
% measured force for each image.

% Normally, for this method, we have to launch many finite element
% computations. However, to speed up this workshop, and to bypass some
% licensing issues. A cloud of possible solutions is precomputed using
% Abaqus. Today, we will only use this data and interpolate our specific
% results from it using the abaqus.m function.

% In FEMU, like in DIC we start by defining a cost function. For example,
% for the distance between the FEM and the experimental displacement
% fields:
% Cu = sum (Ufem - Uexp)^2;
% However, in this example we also have the force, which leads to a second
% shape function
% Cf = sum (Ffem - Fexp)^2;
% The total cost function is the sum of both, but before we can sum them we
% have to scale them accordingly. C = Wu * Cu + Wf * Cf
% We will discuss this scaling.

% Again, like DIC, we apply Newton's method to the cost function to get the
% linearised system of equations that we solve iteratively, starting from
% an initial guess. However, this time the Degrees of Freedom (DOF) are the
% material parameters.

% The material chosen to represent this experiment is linear elasticity
% with an exponential (Voce) hardening law using 5 parameters:
% p(1) : Young modulus
% p(2) : Poisson ratio
% p(3) : Yield Stress
% p(4) : Ultimate Stress
% p(5) : Hardening coefficient

% the constitutive law writes:
% sigma = p(1) * epsilon, if sigma < p(3)
% sigma = p(3) + (p(4) - p(3)) * ( 1 - exp( -p(5) * epsilon_p ));

% Let's start programming
% --------------------------------------

% some strings giving names to our 5 parameters
matstr = {'Young modulus','Poisson ratio','Yield Stress',...
 'Ultimate Stress','Hardening coefficient'};

% an initial for the five parameters
p0 = [70000; 0.3; 400; 650; 70];
p = p0;

% load the precomputed data
ABQ = load('ABQ.mat');

% get the pixel size to relate pixels to [mm]
pxsize = ABQ.pxsize;

% get the mesh
coor = %YOUR CODE HERE
conn = %YOUR CODE HERE

% get the experimental data in real units
Uexp = %YOUR CODE HERE
Fexp = %YOUR CODE HERE

% get some matrix sizes
Nn = size(coor,1);
Ne = size(conn,1);
Ninc = size(Fexp,2);
Np = numel(p);

% indices to the separate displacement components
Ix = 1:2:2*Nn;
Iy = 2:2:2*Nn;

% an estimate for the measurement uncertainty
gamma_u = 0.01; % px
gamma_f = 0.1; % N

% Prepare a figure for plotting
figure('Position',[50 50 1400 800]);
for inc = 1:Ninc
    ha(inc) = subplot(2,Ninc,inc);
    hp(inc) = patch('Vertices',coor,'Faces',conn,'FaceVertexCData',zeros(Nn,1),'EdgeColor','k','FaceColor','interp');
    title(sprintf('Ru [mm] - inc(%d)',inc));
    colorbar
end
ha(inc+1) = subplot(2,Ninc,Ninc + (1:Ninc),'NextPlot','add');
hp(inc+1) = plot(1:Ninc,zeros(Ninc,1));
hp(inc+2) = plot(1:Ninc,zeros(Ninc,1),'Color',0.5 * [1, 1, 1]);
title(sprintf('Rf',inc));
xlabel('increments');
ylabel('rF [N]');

% maximum number of iterations
maxit = 20;

% convergence criteria (on dp)
convcrit = 1e-3;

history.p = zeros(maxit,Np);
history.R = zeros(maxit,3);
history.Rf = zeros(maxit,Ninc);

% printing a header line
fprintf('   %3s, %10s, %10s, %10s, %10s\n','it','R','Rf','Ru','dp');

% the FEMU loop
for it = 1:maxit

    % compute the current solution and the sensitivity fields
    [Uref, Fref, Su, Sf] = abaqus(p, ABQ);

    % compute the residual
    Ru = %YOUR CODE HERE
    Rf = %YOUR CODE HERE

    % the total residual
    R = sqrt(Wu) * rms(Ru) + sqrt(Wf) * rms(Rf);

    % compute the Hessian
    Hu = %YOUR CODE HERE
    Hf = %YOUR CODE HERE

    % compute the Jacobian matrix
    Ju = %YOUR CODE HERE
    Jf = %YOUR CODE HERE

    % combine the two systems
    H = %YOUR CODE HERE
    J = %YOUR CODE HERE

    % a little bit of Levenbergh-Marquardt
    alpha = 0.001 * eigs(H,1);
    H = H + diag(alpha) .* eye(Np);

    % solve for the update in the parameters:
    dp = %YOUR CODE HERE

    % update the parameters
    p = %YOUR CODE HERE

    % printing progress
    fprintf('   %3d, %10.3e, %10.3e, %10.3e, %10.3e\n',it,R,rms(Rf),rms(Ru),rms(dp));

    % update some history
    history.p(it,:) = p;
    history.R(it,:) = [R,rms(Rf),rms(Ru)];
    history.Rf(it,:) = Rf;

    % update the figure
    for inc = 1:Ninc
        % the offset
        O = (inc - 1) * 2 * Nn;

        % the residual amplitude
        Ra = sqrt( Ru(O + Ix).^2 + Ru(O + Iy).^2 );
        set(hp(inc),'FaceVertexCData',Ra);
    end
    set(hp(inc+1),'YData',Rf);
    set(hp(inc+2),'YData',history.Rf(1,:));
    drawnow
end

savepng('residual_FEMU.png')

% Print the found parameters to the screen
fprintf('The identified parameters:\n')
for k = 1:Np
    fprintf('%s : %10.3e\n',matstr{k}, p(k));
end

% Show the parameter convergence behaviour
figure;
plot(1:maxit, history.p ./ transpose(p0), '.-' );
legend(matstr)

savepng('convergence_FEMU.png')