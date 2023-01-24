clear; close all

% Exercise 1 : DIC
% -----------------------------------------------------------
% DIC is an easy starting point for the discussion about inverse methods.
% However, it is not a necessary one. DIC is just an example of an inverse
% method.

% Nevertheless, for the sake of this workshop we start with DIC as the form
% a nice warm-up exercise getting familiar with many of the relevant
% methods and equations.

% There are many flavours of DIC, each with their own theoretical advantages
% that if used properly are less significant in practice. The version used
% today is a Finite Element (FE) based, Global DIC (GDIC), code that uses
% T3 elements. The main advantage of this method is that it is easy to
% connect it to a FE analyses since they use the same basis.

% load some functions provided in the lib folder
addpath('lib');

% Images
% -----------------------------------------------------------

% Define the path to the images
imagepath = 'images';

% getting a list of all images
imagefiles = dir(fullfile(imagepath, '*.tif'));
imagefiles = {imagefiles.name}';
images = strcat(imagepath, filesep, imagefiles);

% the number of increments (experimental steps)
Ninc = numel(images) - 1;

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

% !! why convert to float?

% Mesh
% -----------------------------------------------------------

load('DIC_01.mat')

Nn = size(coor,1);
Ne = size(conn,1);

figure;
imagesc(frgb);
patch('Vertices',coor,'Faces',conn,'EdgeColor','y','FaceColor','none','LineWidth',2);
axis image ij
colorbar

% DIC - pair of images
% -----------------------------------------------------------

% Let's write a DIC code from scratch, it is not that hard. We want to
% minimise the brightness conservation equation: f(x) - g(x + u(x)) = 0

% With a bit of math, we end up with a linearised system of equations that
% we will iteratively solve to approach towards to solution using Newton's
% method. The solution are the nodal displacement amplitudes stored as "a",
% which after initialising are the iteratively updated with:
% a = a + da

% To find the iterative update "da" we have to solve the system:
% M * da = b --> da = inv(M) * B, where
% M_ij = sum (gradf * phi_i) * (phi_j * gradf), and
% b_i = sum r * (phi_i * gradf)

% notice that this system only has three ingredients:
% - gradf -> the image gradient
% - phi -> the shapefunction matrix
% - r -> the residual image: r = (f(x) - g(x + u(x)))

% Now in code
% -----------------

% compute the shapefunction matrix
[phi, E] = TriangleShapefunGridded(coor, conn, size(f));

% find all pixels inside the mesh
mask = E > 0;

% show one shapefunction
shownode = 49;
figure;
imagesc( reshape(phi(:,shownode), n, m) );
patch('Vertices',coor,'Faces',conn,'EdgeColor','k','FaceColor','none','LineWidth',2);
axis image ij
colorbar

% compute the gradient images
[fx, fy] = gradient(f);

% show the gradient images
figure;
subplot(1,2,1)
colormap(gray);
imagesc(fx);
title('fx')
axis image ij
colorbar

subplot(1,2,2)
colormap(gray);
imagesc(fy);
title('fy')
axis image ij
colorbar

% DOF indices (the x and y degrees of freedom are interleaved)
Ix = 1:2:2*Nn;
Iy = 2:2:2*Nn;

% Hessian (the transpose multiplication is performing the sum)
M = %YOUR CODE HERE

% show the hessian
figure;
spy(M)

% load the second image
g = %YOUR CODE HERE

% prepare the interpolator
interpmeth = 'linear';
extrapval = nan;

% initialize the displacement field
a = %YOUR CODE HERE

% show the initial residual
figure;
imagesc(f - g)
patch('Vertices',coor,'Faces',conn,'EdgeColor','y','FaceColor','none');
colormap(gray)
title('initial residual')
axis image ij
caxis(5e3*[-1, 1]);
colorbar

% maximum number of iterations
maxit = 20;

% convergence criteria (on da)
convcrit = 1e-3;

% printing a header line
fprintf('   %3s, %10s, %10s\n','it','res','da');

% start the DIC loop
for it = 1:maxit

    % get the current displacement for each pixel
    Ux = %YOUR CODE HERE
    Uy = %YOUR CODE HERE

    % interpolate the image at the deformed locations (g-tilde)
    gt = %YOUR CODE HERE

    % compute the residual image
    res = %YOUR CODE HERE

    % create the Jacobian matrix
    b = %YOUR CODE HERE

    % solve the system
    da = %YOUR CODE HERE

    % update the dof
    a = %YOUR CODE HERE

    % compute the rms of this
    r = rms(res(mask));

    % print some status
    fprintf('   %3d, %10.3e, %10.3e\n',it,r,rms(da));

    if rms(da) < convcrit
        fprintf('converged\n\n')
        break
    end
end

% show the final residual
figure;
imagesc(res)
patch('Vertices',coor,'Faces',conn,'EdgeColor','y','FaceColor','none');
colormap(gray)
title('final residual')
axis image ij
caxis(5e3*[-1, 1]);
colorbar


% DIC - series of images
% -----------------------------------------------------------

% The function dic_T3 encapsulates all the above code in a single function
% such that we can compactly write a code that runs on a series of images.
% It also implements this idea of a pyramid approach where DIC is first
% performed on blurred images, the solution of which is used as initial
% guess for less blurred images etc. Until the final pyramid step that
% performs DIC on the untouched images.

% options
opt.blur = [10, 5, 1, 0];

% initial guess
init = %YOUR CODE HERE

for inc = 1:Ninc

    fprintf('Increment %d/%d\n',inc,Ninc);
    fprintf('======================================\n');


    % read the current image
    g = %YOUR CODE HERE

    % perform DIC
    cor(inc) = dic_T3(f,g,init,coor,conn,opt);

    % update the initial guess
    init = %YOUR CODE HERE

    % re-use the figure
    opt.Hax = gca;
end
