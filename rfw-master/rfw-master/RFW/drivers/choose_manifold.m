function fs = choose_manifold(manifold, retraction, vector_transport)
% This function return some necessary function handles of some manifolds which are required in Riemannian optimization methods.
%
% INPUT:
% retraction : a number represents a retraction. See sphere.m, stiefel.m orthogroup.m and grassmann.m for details
% vector transport : a number represents a vector transport. See sphere.m, stiefel.m orthogroup.m and grassmann.m for details
% OUTPUT:
% Output function handles about the SPD matrix manifold.

fs = SPD(retraction, vector_transport);

end
