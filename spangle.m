function C=spangle_3(u,V)
% u are vectors scaled 0-1 and should be same class.
% C, V are matrices scaled same.

% u=u(:);
%% From Mathworks answers
[a,b,c]=size(V);
u=permute(u, [3 1 2]);
U=repmat(u, [a, b, 1]);
Cr = cross(U,V);
NoCr= sqrt(sum(Cr.*Cr, 3));%norm of cross
C = atan2d(NoCr, dot(U,V, 3));


