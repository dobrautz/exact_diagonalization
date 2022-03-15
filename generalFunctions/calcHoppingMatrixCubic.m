function hoppingMatrix = calcHoppingMatrixCubic(lattice, param)
%calcHoppingMatrix calculates the hopping matrix for a 2 dimensional cubic
% lattice for periodic and twisted boundary conditions
% 
% Input:
% 
% Output:
% 
% Notes: as title says for now only correct for a dimensional cubic lattice

% hopping matrix
hoppingMatrix = zeros(lattice.nSites);
% set connected sites to -t
linIndexNeighbors = sub2ind(size(hoppingMatrix),repmat(1:lattice.nSites,lattice.dim,1),...
    lattice.neighbors);
linIndNeighRight = linIndexNeighbors(1,:);
linIndNeighDown = linIndexNeighbors(2,:);

hoppingMatrix(linIndNeighRight) = -param.t * exp(1i*lattice.twisted(1));
% negaitive phase for down hopping
hoppingMatrix(linIndNeighDown) = -param.t * exp(-1i*lattice.twisted(2));

% since t is symmetric but neighors are only in left and down direction:
hoppingMatrix = hoppingMatrix + hoppingMatrix';