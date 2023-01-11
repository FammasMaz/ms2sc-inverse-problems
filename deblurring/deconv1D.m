% Blurring operator (G. Stadler, U. of Texas)

% Number of discretization points
N = 64;
% N = 128;
gamma = 0.03;
C = 1 / (sqrt(2*pi)*gamma);

K = zeros(N,N);
h = 1/N;
x = linspace(0,1,N)';

% discrete convolution matrix
K = h * C * exp(-(toeplitz(0:N-1)).^2 * h^2 / (2 * gamma^2));
% alternative way to calculate K
%for l = 1:N
%    for k = 1:N
%    	K(l,k) = h * C * exp(-(l-k)^2 * h^2 / (2 * gamma^2));
%    end
%end

% exact parameters
p = (x > .2).*(x < .3) + sin(4*pi*x).*(x > 0.5);

% convolved parameters
d = K * p;

% noisy data, noise has sigma^2 = 0.01
dn = d + 0.1 * randn(N,1);
figure; plot(x,d,x,dn,x,p,'Linewidth', 2);
legend('data', 'noisy data', 'exact parameters');

% direct inverse with perfect data
p1 = inv(K) * d;
figure; plot(x,p1,x,p,'Linewidth', 2);
legend('parameters with perfect data', 'exact parameters');

% direct inverse with noisy data
p2 = inv(K) * dn ;
figure; plot(x,p2,x,p,'Linewidth', 2);
legend('parameters with noisy data', 'exact parameters');
