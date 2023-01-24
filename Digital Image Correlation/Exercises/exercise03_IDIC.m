clear; close all

% load some functions provided in the lib folder
addpath('lib');

% Exercise 3 : Weighted FEMU
% -----------------------------------------------------------

% In the previous exercise, we gave an equal weight to all the nodes. This
% is wrong. The elements were of different sizes. Consequently, the support
% of each node (the area of all its connected elements) was not the same.
% Moreover, not every part of the image contained in equal amount of
% information.

% The DIC "M-matrix" contains all the information that is required to
% properly weight the contribution to each node. This has been called
% M-weighted FEMU.

% Integrated DIC (IDIC) is a method that combines FEMU and DIC into a
% single step. This naturally includes the M-weighting. As was shown in the
% lecture, M-weighted FEMU and IDIC are mathematically equivalent. However,
% the main advantage of IDIC is that it strongly regularises the DIC
% method. Consequently, we can adopt a much finer mesh and thus have more
% accurate results in the DIC part.


% A repeat of the FEMU code from exercise 2
% -------------------------------------------------------

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

% an estimate for the measurement uncertainty
gamma_u = 0.01; % px
gamma_f = 0.1; % N
gamma_i = 1000; % grayvalues

% Weighted FEMU - using code from exercise 1
% -------------------------------------------------------

% Lets get the DIC M-matrix

% Define the path to the images
imagepath = 'images';

% getting a list of all images
imagefiles = dir(fullfile(imagepath, '*.tif'));
imagefiles = {imagefiles.name}';
images = strcat(imagepath, filesep, imagefiles);

% Load the first image
f = imread(images{1});

% create an RGB version of f for later plotting
frgb = repmat(f, 1, 1, 3);

dynamicrange = max(f(:)) - min(f(:));
fprintf('Dynamic range of f: %g [GV] \n', dynamicrange);

% get the size of the image
[n, m] = size(f);

% convert to floating point
f = double(f);

% compute the gradient images
[fx, fy] = gradient(f);

% compute the shapefunction matrix
[phi, E] = TriangleShapefunGridded(ABQ.coor, conn, size(f));

% find all pixels inside the mesh
mask = E > 0;

% Hessian (the transpose multiplication is performing the sum)
M = %YOUR CODE HERE

% Integrated DIC - using code from exercise 1
% -------------------------------------------------------

% prepare the interpolator
interpmeth = 'linear';
extrapval = nan;


% The loop
% -------------------------------------------------------

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
fprintf('   %3s, %10s, %10s, %10s, %12s, %10s\n','it','R','Rf','Ru','Ri','dp');

% switch to use FEMU or IDIC
use_FEMU = true;

% the FEMU / IDIC loop
for it = 1:maxit

    % compute the current solution and the sensitivity fields
    [Uref, Fref, Su, Sf] = abaqus(p, ABQ);
    
    % compute the residual
    Ru = %YOUR CODE HERE
    Rf = %YOUR CODE HERE
    
    % the total residual
    R = sqrt(Wu) * rms(Ru) + sqrt(Wf) * rms(Rf);

    % introduce the DIC weighting per increment for the displacements
    if use_FEMU
        Hu = %YOUR CODE HERE
        Ju = %YOUR CODE HERE
    else % IDIC
        Hu = %YOUR CODE HERE
        Ju = %YOUR CODE HERE        
    end
    
    % the force Hessian and Jacobian matrix
    Hf = %YOUR CODE HERE
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
    fprintf('   %3d, %10.3e, %10.3e, %10.3e, %12.5e, %10.3e\n',it,R,rms(Rf),rms(Ru),rms(Ri),rms(dp));

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

if use_FEMU
    savepng('residual_FEMU_M.png')
else
    savepng('residual_IDIC.png')
end


% Print the found parameters to the screen
fprintf('The identified parameters:\n')
for k = 1:Np
    fprintf('%s : %10.3e\n',matstr{k}, p(k));
end

% Show the parameter convergence behaviour
figure;
plot(1:maxit, history.p ./ transpose(p0), '.-' );
legend(matstr)


if use_FEMU
    savepng('convergence_FEMU_M.png')
else
    savepng('convergence_IDIC.png')
end




