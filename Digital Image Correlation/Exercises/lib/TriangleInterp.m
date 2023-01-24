function a2 = TriangleInterp(coor1,a1,coor2)
% a2 = TriangleInterp(coor1,a1,coor2), convert the output
% obtained using one mesh (coor1) to the coordinates of a second mesh
% as defined by coor2.
%
% This function will perform a natural interpolation and linar
% extrapolation using the scatteredInterpolant matlab function

Nn = size(coor2,1);
a2 = zeros(Nn,2);

G = scatteredInterpolant(coor1,a1(:,1),'natural','linear');
a2(:,1) = G(coor2);
G.Values = a1(:,2);
a2(:,2) = G(coor2);



