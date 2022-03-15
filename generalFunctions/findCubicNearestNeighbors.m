function [ indNeighbors , notPeriodicBound] = findCubicNearestNeighbors( latticeDim )
%findCubicNearestNeighbors gives indices of nearest neighbors for cubic
%   lattices in 1, 2 and 3 Dimensions
% 
%   indNeighbors  = findCubicNearestNeighbors( latticeDim ) assumes 
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
% 
%------------------------SVN Info------------------------------------------
% $Rev:: 51                             $: Revision of last commit
% $Author:: dobrautz                    $: Author of last commit 
% $Date:: 2013-05-30 21:28:42 +0200 (Do#$: Date of last commit
% -------------------------------------------------------------------------

% periodic boundary conditions are assumed to be standard
[absoluteDim,~] = size(latticeDim);   %dimension des problems
nSites = prod(latticeDim , 1);          %Laenge der Heisenbergkette

% control matric for function circshift() for different dimensions
circshiftControl = [0,-1,0;-1,0,0;0,0,-1];     

indNeighbors = zeros([absoluteDim,latticeDim']);

if absoluteDim > 1  % check for dimensions
    index = reshape(1:nSites , latticeDim');
else
    index = 1:nSites;
end

for l = 1:absoluteDim % loop over dimensions
    indNeighbors(l,:,:,:) = circshift(index , circshiftControl(l,:));
end

indNeighbors = reshape(indNeighbors,[absoluteDim,numel(indNeighbors)/absoluteDim]);

% nearest neighbours nicht periodische RB
if nargout > 1 
    notPeriodicBound = zeros([absoluteDim,latticeDim']);
    if absoluteDim == 1
        index = [1:nSites,0];
        
        for l = 1:absoluteDim
            a = circshift(index,circshiftControl(l,:));
            notPeriodicBound(l,:,:,:) = a(1:latticeDim(l));
        end
        
    elseif absoluteDim ==2
        index = cat(1,reshape(1:nSites,latticeDim'),zeros(1,latticeDim(2)));
        index = cat(2,index,zeros(latticeDim(1)+1,1));
        for l = 1:absoluteDim
            a = circshift(index,circshiftControl(l,:));
            notPeriodicBound(l,:,:) = a(1:latticeDim(1),1:latticeDim(2));
        end
        
    elseif absoluteDim == 3
        index = cat(1,reshape(1:nSites,latticeDim'),zeros(1,latticeDim(2),latticeDim(3)));
        index = cat(2,index,zeros(latticeDim(1)+1,1,latticeDim(3)));
        index = cat(3,index,zeros(latticeDim(1)+1,latticeDim(2)+1,1));
        
        for l = 1:absoluteDim
            a = circshift(index,circshiftControl(l,:));
            notPeriodicBound(l,:,:,:) = a(1:latticeDim(1),1:latticeDim(2),1:latticeDim(3));
        end
    end
    
    notPeriodicBound = reshape(notPeriodicBound,[absoluteDim,numel(notPeriodicBound)/absoluteDim]);
    
end

end

