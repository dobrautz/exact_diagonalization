function [ indNeighbors , nSites, notPeriodicBound ] = getNearestNeighbors(varargin)
%GETNEARESTNEIGHBORS calculates the nearest neighbor matrix for different
%standard lattices: cubic, square-root cubic, triangular, honeycomb and
%kagome lattice.
% Input:        latticeType         string choosing the lattice type
%               latticeSize         size of desired lattice
%
% Output:       indNeighbors        index matrix of nearest neighbors
%               nSites              number of Sites
%               flagNonPeriodic     flag indicating non oeriodic boundary
%                                   conditions
%
% NOTE: 13.07.13: first draft of this function
%
%------------------------SVN Info------------------------------------------
% $Rev::                                        $: Revision of last commit
% $Author::                                     $: Author of last commit
% $Date::                                       $: Date of last commit
% -------------------------------------------------------------------------

if numel(varargin) == 1

    latticeSize = varargin{1}.size;
    latticeType = varargin{1}.type;
    
elseif numel(varargin) == 2
    
    latticeType = varargin{1};
    latticeSize = varargin{2};
    
end


[latticeDim, ~] = size(latticeSize);

switch latticeType
    
    case 'cubic'
        
        nSites = prod(latticeSize , 1);          %Laenge der Heisenbergkette
        
        % control matric for function circshift() for different dimensions
        circshiftControl = [0,-1,0;-1,0,0;0,0,-1];
        
        indNeighbors = zeros([latticeDim,latticeSize']);
        
        if latticeDim > 1  % check for dimensions
            index = reshape(1:nSites , latticeSize');
        else
            index = 1:nSites;
        end
        
        for l = 1:latticeDim % loop over dimensions
            indNeighbors(l,:,:,:) = circshift(index , circshiftControl(l,:));
        end
        
        indNeighbors = reshape(indNeighbors,[latticeDim,numel(indNeighbors)/latticeDim]);
        
        % nearest neighbours nicht periodische RB
        if nargout > 1
            notPeriodicBound = zeros([latticeDim,latticeSize']);
            if latticeDim == 1
                index = [1:nSites,0];
                
                for l = 1:latticeDim
                    a = circshift(index,circshiftControl(l,:));
                    notPeriodicBound(l,:,:,:) = a(1:latticeSize(l));
                end
                
            elseif latticeDim ==2
                index = cat(1,reshape(1:nSites,latticeSize'),zeros(1,latticeSize(2)));
                index = cat(2,index,zeros(latticeSize(1)+1,1));
                for l = 1:latticeDim
                    a = circshift(index,circshiftControl(l,:));
                    notPeriodicBound(l,:,:) = a(1:latticeSize(1),1:latticeSize(2));
                end
                
            elseif latticeDim == 3
                index = cat(1,reshape(1:nSites,latticeSize'),zeros(1,latticeSize(2),latticeSize(3)));
                index = cat(2,index,zeros(latticeSize(1)+1,1,latticeSize(3)));
                index = cat(3,index,zeros(latticeSize(1)+1,latticeSize(2)+1,1));
                
                for l = 1:latticeDim
                    a = circshift(index,circshiftControl(l,:));
                    notPeriodicBound(l,:,:,:) = a(1:latticeSize(1),1:latticeSize(2),1:latticeSize(3));
                end
            end
            
            notPeriodicBound = reshape(notPeriodicBound,[latticeDim,numel(notPeriodicBound)/latticeDim]);
            
        end
        
        
        
    case 'squareRoot'
        
        %TODO: allgemeiner!
        % number of Sites:2D eg. = Lx*Ly + (Lx-1)*(Ly-1)
        nSites = prod(latticeSize,1) + prod(latticeSize-1,1);
        
        % first idea how get nearest neigbor relation: index states in such a way
        % so nearest neighbor are easy not other way around
        
        indNeighbors = zeros(latticeDim,nSites);
        
        siteIndex = (1:nSites);
        % right neighbor always + 1 to current index
        indNeighbors(1,:) = circshift(siteIndex,[0,-1]);
        
        for iDim = 2:latticeDim
            
            if sum(latticeSize) ~= 7
                indNeighbors(iDim,:) = circshift(siteIndex,[0,(min(latticeSize)-1)^2+1]);
                
            else
                
                indNeighbors(iDim,:) = circshift(siteIndex,[0,8]);
                
                indNeighbors(1,:) = circshift(siteIndex,[0,-3]);
                
            end
            
        end
    case 'triangular'
        % NOTE: 13.07.13: for now only 'quadratic' triangular! and only 2D
        % and just PBC
        
        nSites = prod(latticeSize);
        indSites = 1:nSites;
        
        % latticeDim+1 because 6 nearest neigbors in 2D triangular
        indNeighbors = zeros(latticeDim+1,nSites);
        % first 2 rows of NN matrix will be translaitonal directions(down-right)
        % and to the right) = same as in quadratic cubic lattice
        circshiftControl = [0,-1,0;-1,0,0;1,-1,0];
                
        indSites2D = reshape(indSites , latticeSize');
          
        % first 2 NN same as quadratic
        % 3rd is a shift to the left and down!
        for iDim = 1:latticeDim+1 % loop over dimensions
            transNN = circshift(indSites2D , circshiftControl(iDim,:));
            
            indNeighbors(iDim,:) = transNN(:);
            
        end
        
    case 'honeycomb'
        
        %TODO
    case 'kagome'
        
        %TODO
    case 'tilted'
        
        % double so many lattice points:
        nSites = 2*prod(latticeSize);
        
        % neighbors probably similar to square root one TODO!!!
%         NOTE: 28.10.13: found solution see documentation
        siteIndex = 1:(2*latticeSize(1)*latticeSize(2));
        siteIndex = reshape(siteIndex,latticeSize(1),2*latticeSize(2));
        
        indNNright = circshift(siteIndex,[0,-1]);
        
        indNNright(:,2:2:end) = circshift(indNNright(:,2:2:end),-1);
        
        indNNdown = circshift(indNNright(:),2*latticeSize(1));
        
        indNeighbors = [indNNright(:)';indNNdown'];
        
end

