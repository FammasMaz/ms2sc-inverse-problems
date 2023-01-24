function [U, F, Su, Sf] = abaqus(p, ABQ)
% [U, F] = abaqus(p, ABQ), Compute the displacement field and reaction
% force for a given set of parameters p using the precomputed database ABQ
%
% [U, F, SU, SF] = abaqus(p, ABQ), Also return it's derivatives, i.e. the
% sensitivity fields. The returned sensitivity fields are computed to the
% relative parameters i.e. (Uper - Uref) / perfactor, instead of the
% absolute version (Uper - Uref) / (p_per - p_ref).

T = ABQ.T;
P = ABQ.P;

Nt = size(T,1);
Np = size(P,1);

d = size(P,2);
n = size(T,2);

% Compute S?
wantS = nargout > 2;

% force p to be a row vector
p = transpose(p(:));

% normalize P
% --------------------------------------

Pn = P ./ p;
pn = p ./ p;

% get the centers of the elements
% --------------------------------------
% this is not required, but it speads up the search in the loop below

% element centers
T1 = Pn(T(:,1));
T2 = Pn(T(:,2));
T3 = Pn(T(:,3));
T4 = Pn(T(:,4));
T5 = Pn(T(:,5));
T6 = Pn(T(:,6));
Pc = (1/6) * (T1 + T2 + T3 + T4 + T5 + T6);

% compute the distnace from p to each element
dP = sum( (Pc - pn).^2, 2);

% sort them in order
[~, searchorder] = sort(dP);


% test each element to see of p is inside
% --------------------------------------

mar = 1e-5;
best.t = 0;
best.r = Inf;
found_t = false;
for kt = 1:Nt
    
    % get the element number
    t = searchorder(kt);
    
    % get the one element
    con = T(t,:);
    
    % get the element coordinates
    x = Pn(con,:);

    % element volume
    V = det([ones(n,1),x]);
    
    % compute the shapefunctions
    M = cat(2, x(1:d+1,1:d), ones(d+1,1));
    b = cat(2, pn, 1);
    L = b / M;    
    
    % test if p is inside t
    if all(L > (0 - mar) & L < (1+mar))
        found_t = true;
        % found the correct simplex
        break
    else
        % compute the distance from the center
        r = sum( (L - 0.5).^2 );
        if r < best.r
            best.r = r;
            best.t = t;
        end
    end
end

% if we didn't find an element, pick the nearest one
if ~found_t
    t = best.t;
end
assert(t ~= 0,'solution not found, parameters are out of the reasonable space');

% Interpolate the results
% --------------------------------------

% get the one element
con = T(t,:);

% get the derivatives
if wantS
    % get the element coordinates
    x = Pn(con,:);
    
    % Jacobian matrix
    J = x(1:end-1,:) - x(end,:);
    
    % Gradient
    dL = inv(transpose(J));
    
    % add the missing grad -> L(6) = 1 - sum(L(1:5))
    dL(end+1,:) = -sum(dL,1);
end

% prepare storage
Nu = numel(ABQ.Uexp);
Nf = numel(ABQ.Fexp);
U = zeros(Nu,1);
F = zeros(Nf,1);
if wantS
    Su = zeros(Nu,5);
    Sf = zeros(Nf,5);
end

% interpolate
for k = 1:n
    
    % get the data for this vertex (as double)
    Uk = double(ABQ.fem(con(k)).U(:));
    Fk = double(ABQ.fem(con(k)).F(:));
    
    % interpolate
    U = U + L(k) .* Uk;
    F = F + L(k) .* Fk;
    
    % also get it's derivative
    if wantS
        Su = Su + dL(k,:) .* Uk;
        Sf = Sf + dL(k,:) .* Fk;
    end
end
