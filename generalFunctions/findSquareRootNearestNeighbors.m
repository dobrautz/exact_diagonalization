function [ indNeighbors , notPeriodicBound] = findSquareRootNearestNeighbors(latticeDim)
%FINDSQUAREROOTNEARESTNEIGHBORS gives indices of nearest neighbors for 
% 'square-root' lattice in 2 and 3 dimensions
% 
% indNeighbors  = findSquareRootNearestNeighbors( latticeDim ) assumes 
%   periodic boundary conditions
% 
%   [indNeighbors, notPeriodicBound] = findCubicNearestNeighbors(latticeDim)
%   a call with two output arguments assumes non periodic boundary conditions
%   for which a flag index output notPeriodicBound is also given
%
%   Input:      latticeDim   ...     lattice dimensions
%               
%   Output:     indNeighbors     ... Indizes nearest Neighbours bei 
%                                    periodischen RB
%               notPeriodicBound ... Kombination mit nn Indizes bei nicht 
%                                    per. RB
% NOTE: 10.07.13: first draft of this function: 
%  the number of lattice sites is determined by formula:
% nSites = Lx*Ly + (Lx-1)*(Ly-1) easy to verify! where Lx, Ly is the number
% of next nearest neighbors on the edges of the lattice
% Indexing of sites in such manner:     2
%                                     1 3 5
%                                       4  NO! arrange indizes so nn
%                                        indexing is easy! see below
% NOTE: 10.07.13: many much TODO
%------------------------SVN Info------------------------------------------
% $Rev::                                        $: Revision of last commit
% $Author::                                     $: Author of last commit 
% $Date::                                       $: Date of last commit
% -------------------------------------------------------------------------

% keyboard
% periodic boundary conditions are assumed to be standard
[absoluteDim,~] = size(latticeDim);   %dimension des problems
% number of Sites:2D eg. = Lx*Ly + (Lx-1)*(Ly-1)
nSites = prod(latticeDim,1) + prod(latticeDim-1,1);         

% first idea how get nearest neigbor relation: index states in such a way
% so nearest neighbor are easy not other way around

indNeighbors = zeros(absoluteDim,nSites);

siteIndex = (1:nSites);
% right neighbor always + 1 to current index
indNeighbors(1,:) = circshift(siteIndex,[0,-1]);

for iDim = 2:absoluteDim
    
    if sum(latticeDim) ~= 7
        indNeighbors(iDim,:) = circshift(siteIndex,[0,(min(latticeDim)-1)^2+1]);
    
    else
        
        indNeighbors(iDim,:) = circshift(siteIndex,[0,8]);
        
        indNeighbors(1,:) = circshift(siteIndex,[0,-3]);
        
    end
    
end


unique(indNeighbors(2,:))

test = zeros(10,nSites);
test(1,:) = siteIndex;

for i = 2:10

test(i,:) = indNeighbors(2,test(i-1,:));

end

test