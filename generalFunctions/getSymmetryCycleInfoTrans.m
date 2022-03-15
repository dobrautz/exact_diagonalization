function [ phasePrefactors, indOccupiedUp, indOccupiedDown] = ...
    getSymmetryCycleInfoTrans( basisStateUp, basisStateDown, latticeSize, kValue)
%GETSYMM determines the phases(translation, symmetry,
% fermionic) of members of a cycle for a specific representative for only
% translational symmetry
%
% Input:    basisStatesUp/Down      up and down spin part of basis state
%           symmetry                to be applied symmetries
%
% Output:   phasePrefactor          combined phase prefactors (translation,
%                                   symmetry, fermionic)
%           indOccupiedUp/Down      index of occupied sites
%
% NOTE: 03.08.13: created

% size of Basis and number of Sites (number of sites = number of trans. syms!)
nSites = length(basisStateUp);
% Dimension of lattice
[latticeDim,~] = size(latticeSize);
% 1 dimensional lattice sizes: change to rect. 25.06.13:
nSitesY = latticeSize(1);
% index for occupied sites identification
indSites = 1:nSites;
% number of occupied sites:
nOccupiedUp = sum(basisStateUp);
nOccupiedDown = sum(basisStateDown);
% container variable for new states(symmetry applied to them)
symmetricBasisUp = zeros(nSites,nSites);
symmetricBasisDown= zeros(nSites,nSites);
% cont. var. for respective phase factor due to symmetry exchanges(fermions)
symmetryPhases = ones(nSites,1);

%% all translations on basis
if latticeDim == 1
    
    % calc. of exp phase factor for !careful for sign in exponential!
    exponentialPhases = exp(1i*kValue*(0:nSites-1))';
    
    for rx1D = 0:nSites-1
        % application of Translational symmetry T :
        % index transformation for diff R values
        ind = rx1D + 1;
        
        % application on TransSymx to UP and DOWN spin part:
        [symmetricBasisUp(ind,:), symmetryPhasesUp ] = ...
            cubicSymmetries(basisStateUp, rx1D, latticeSize);
        [symmetricBasisDown(ind,:), symmetryPhasesDown ] = ...
            cubicSymmetries(basisStateDown, rx1D, latticeSize );
        % combine fermionic phases
        symmetryPhases(ind) = symmetryPhasesUp.*symmetryPhasesDown;
        
    end
    
elseif latticeDim == 2
    
    nSitesX = latticeSize(2);
    
    % exp. phases pre-determined from lattice size and symmetry size:
    %-----------------------exponential phases---------------------------------
    rVectorX = 0:nSitesX-1;
    rVectorY = 0:nSitesY-1;
    [rVectorY, rVectorX] = ndgrid(rVectorY, rVectorX);
    rVectorX = rVectorX(:);
    rVectorY = rVectorY(:);
    exponentialPhases = exp(1i*(kValue(1)*rVectorX + kValue(2)*rVectorY));
    % surpress real and imaginary parts of this below treshhold
    exponentialPhases= real(exponentialPhases).*(abs(real(exponentialPhases))>10^-5) + ...
        1i*imag(exponentialPhases).*(abs(imag(exponentialPhases))>10^-5);
    % repmat it to whole symmetry basis ( actually not needed but for
    % uniform indication useful
%     exponentialPhases = repmat(exponentialPhases,combSymmetries,1);
    
    
    for rx = 0:nSitesX-1
        for ry = 0:nSitesY-1
            % application of Translational symmetry T :
            % index transformation for diff R values
            ind = nSitesY*rx + ry + 1;
            % UP SPIN PART: application on Translational Symmetry X
            [ symAppliedStates , phasesUpX ] = cubicSymmetries( basisStateUp , rx, latticeSize );
            % Ty
            [ symmetricBasisUp(ind,:) , phasesUpY ] = cubicSymmetries( symAppliedStates , ry, latticeSize );
            % DOWN SPIN PART: application on Translational Symmetry X
            [ symAppliedStates , phasesDownX ] = cubicSymmetries( basisStateDown , rx, latticeSize );
            % Ty
            [ symmetricBasisDown(ind,:) , phasesDownY ] = cubicSymmetries( symAppliedStates , ry, latticeSize );
            % save fermionic phases
            symmetryPhases(ind) = phasesUpX .* phasesUpY .* phasesDownX .* phasesDownY;
            
        end
        
    end
end

% combine fermionic transl and symmetryphases/eigenvalue
phasePrefactors = symmetryPhases .* exponentialPhases;

% get indices of occupied sites for up and down spin and bring in form
indSites = repmat(indSites',1,nSites);

indOccupiedUp = indSites(logical(symmetricBasisUp'));
indOccupiedDown = indSites(logical(symmetricBasisDown'));

indOccupiedUp = reshape(indOccupiedUp, nOccupiedUp, nSites)';
indOccupiedDown = reshape(indOccupiedDown, nOccupiedDown, nSites)';

