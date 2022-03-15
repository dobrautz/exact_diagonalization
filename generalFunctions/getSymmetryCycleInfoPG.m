function [ phasePrefactors, indOccupiedUp, indOccupiedDown ] = ...
    getSymmetryCycleInfoPG(basisStateUp, basisStateDown, whichSymmetry, ...
    kValue, latticeSize, varargin)
%getSymmetryCycleInfo determines the phases(translation, symmetry,
% fermionic) of members of a cycle for a specific representative for point
% group symmetreis
%
% Input:    basisStatesUp/Down      up and down spin part of basis state
%           symmetry                to be applied symmetries
%
% Output:   phasePrefactor          combined phase prefactors (translation,
%                                   symmetry, fermionic)
%           indOccupiedUp/Down      index of occupied sites
%
% NOTE: 03.08.13: created

%% identify which symmetry operation is to be applied
nameFlag_px = strcmp(whichSymmetry,'px');
nameFlag_py = strcmp(whichSymmetry,'py');
nameFlag_pd = strcmp(whichSymmetry,'pd');
nameFlag_pe = strcmp(whichSymmetry,'pe');

% number of variable arguments in = number of symmetry operations!
nSymmetries = numel(varargin);

% check which symmetry has to be applied
if  nSymmetries == 1 && numel(varargin{1}) == 1
    
    if nameFlag_px
        px = varargin{1};
        
    elseif nameFlag_py
        py = varargin{1};
        
    elseif nameFlag_pd
        pd = varargin{1};
        
    elseif nameFlag_pe
        pe = varargin{1};
        
    end
    
elseif nSymmetries == 1 && numel(varargin{1}) == 2
    
    py = varargin{1}(1);
    px = varargin{1}(2);
    
    nSymmetries = 2;
elseif nSymmetries == 1 && numel(varargin{1}) == 3
        
    pd = varargin{1}(1);
    py = varargin{1}(2);
    px = varargin{1}(3);

    nSymmetries = 3;
elseif nSymmetries == 2 % case: two lattice syms
    py = varargin{1};
    px = varargin{2};
    
elseif nSymmetries == 3 % case: 3 lattice syms
    pd = varargin{1};
    py = varargin{2};
    px = varargin{3};

end

%% parameters and container variables
% number of respective combinations PdPy,PdPx... etc.
combSymmetries = 2^nSymmetries;
% number of Sites (1D and 2D) and corr. index
nSites = length(basisStateUp);
indSites = 1:nSites;
nSitesY = latticeSize(1);
nSitesX = latticeSize(2);
% size of symmetry group
sizeSymmetryGroup = combSymmetries*nSites;
% number of occupied sites:
nOccupiedUp = sum(basisStateUp);
nOccupiedDown = sum(basisStateDown);
% container variable for new states(symmetry applioed to them)
symmetricBasisUp = zeros(sizeSymmetryGroup,nSites);
symmetricBasisDown = zeros(sizeSymmetryGroup,nSites);
% cont. var. for respective phase factor due to symmetry exchanges(fermions)
symmetryPhases = ones(sizeSymmetryGroup,1);

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
exponentialPhases = repmat(exponentialPhases,combSymmetries,1);

%--------------------------------------------------------------------------

%% all translations and corresponding point group syms on basis
% loops over all possible R values (square lattice!)
for rx = 0:nSitesX-1
    for ry = 0:nSitesY-1
        
        % application of Translational symmetry T :
        % index transformation for diff R values
        ind = nSitesY*rx + ry + 1;
        % index within Symmetry
        indCycle = ind : sizeSymmetryGroup : (sizeSymmetryGroup - sizeSymmetryGroup + ind);
        
        
        % States get stored in following way within matrix T:
        % T0B1
        % T1B1
        % T2B1
        % ...
        % T0PxB1
        % T1PxB1
        % ...
        % T0PyB1
        % ...
        % T0B2
        % T1B2
        % ...
        
        
        % application of symmetry operation on UP spin part: X first
        [ symAppliedStates , phasesUpX ] = cubicSymmetries( basisStateUp , rx, latticeSize );
        % then Y translation
        [ symmetricBasisUp(indCycle,:) , phasesUpY ] = cubicSymmetries( symAppliedStates , ry, latticeSize );
        % same with DOWN SPIN part:
        [ symAppliedStates , phasesDownX ] = cubicSymmetries( basisStateDown , rx, latticeSize);
        % then Y translation
        [ symmetricBasisDown(indCycle,:) , phasesDownY ] = cubicSymmetries( symAppliedStates , ry, latticeSize);
        
        % save fermionic phases
        symmetryPhases(indCycle) = (phasesUpX .* phasesUpY .* phasesDownX .* phasesDownY) ;
        
        % bis hier immer
        
        % now decide which point group is to be applied:
        if nameFlag_px || nSymmetries > 1
            %application of Translational symmetry T and Px
            indCycle = indCycle + nSites; %increase ind with number of translational syms=x
            
            [ symAppliedStates , ph_p ] = cubicSymmetries( basisStateUp , px, latticeSize);
            
            [ symAppliedStates , phasesUpX ] = cubicSymmetries( symAppliedStates , rx, latticeSize);
            
            [ symmetricBasisUp(indCycle,:) , phasesUpY ] = cubicSymmetries( symAppliedStates , ry, latticeSize);
            
            [ symAppliedStates , ph_p_down ] = cubicSymmetries( basisStateDown , px, latticeSize);
            
            [ symAppliedStates , phasesDownX ] = cubicSymmetries( symAppliedStates , rx, latticeSize);
            
            [ symmetricBasisDown(indCycle,:) , phasesDownY ] = cubicSymmetries( symAppliedStates , ry, latticeSize );
            
            symmetryPhases(indCycle)  = px * (ph_p .* phasesUpX .* phasesUpY .* ...
                ph_p_down .*phasesDownX.*phasesDownY);
            
            
        end
        
        if nameFlag_py || nSymmetries > 1
            %TPy
            indCycle = indCycle + nSites;
            
            [ symAppliedStates , ph_p ] = cubicSymmetries( basisStateUp , py, latticeSize);
            
            [ symAppliedStates , phasesUpX ] = cubicSymmetries( symAppliedStates , rx, latticeSize );
            
            [ symmetricBasisUp(indCycle,:) , phasesUpY ] = cubicSymmetries( symAppliedStates , ry, latticeSize );
            
            [ symAppliedStates , ph_p_down ] = cubicSymmetries( basisStateDown , py, latticeSize);
            
            [ symAppliedStates , phasesDownX ] = cubicSymmetries( symAppliedStates , rx, latticeSize );
            
            [ symmetricBasisDown(indCycle,:) , phasesDownY ] = cubicSymmetries( symAppliedStates , ry, latticeSize );
            
            symmetryPhases(indCycle)  = py * (ph_p .* phasesUpX .* phasesUpY...
                .* ph_p_down .*phasesDownX.*phasesDownY);
            
            
        end
        
        if nSymmetries > 1
            %_________________________________________XY
            
            indCycle = indCycle + nSites;
            
            [ symAppliedStates , ph_px ] = cubicSymmetries( basisStateUp , px, latticeSize);
            
            [ symAppliedStates , ph_py ] = cubicSymmetries( symAppliedStates , py, latticeSize);
            
            [ symAppliedStates , phasesUpX ] = cubicSymmetries( symAppliedStates , rx , latticeSize);
            
            [ symmetricBasisUp(indCycle,:) , phasesUpY ] = cubicSymmetries( symAppliedStates , ry , latticeSize);
            
            [ symAppliedStates , ph_p_down ] = cubicSymmetries( basisStateDown , px, latticeSize);
            
            [ symAppliedStates , ph_py_down ] = cubicSymmetries( symAppliedStates , py, latticeSize);
            
            [ symAppliedStates , phasesDownX ] = cubicSymmetries( symAppliedStates , rx, latticeSize );
            
            [ symmetricBasisDown(indCycle,:) , phasesDownY ] = cubicSymmetries( symAppliedStates , ry, latticeSize );
            
            symmetryPhases(indCycle)  = px*py * (ph_px .* ph_py .* phasesUpX .* ...
                phasesUpY .* ph_p_down .* ph_py_down .* phasesDownX.*phasesDownY);
            
            
        end
        
        if nameFlag_pd || nSymmetries > 2
            %__________________________________D
            
            indCycle = indCycle + nSites;
            
            [ symAppliedStates , ph_p ] = cubicSymmetries( basisStateUp , pd, latticeSize );
            
            [ symAppliedStates , phasesUpX ] = cubicSymmetries( symAppliedStates , rx , latticeSize);
            
            [ symmetricBasisUp(indCycle,:) , phasesUpY ] = cubicSymmetries( symAppliedStates , ry , latticeSize);
            
            [ symAppliedStates , ph_p_down ] = cubicSymmetries( basisStateDown , pd, latticeSize);
            
            [ symAppliedStates , phasesDownX ] = cubicSymmetries( symAppliedStates , rx, latticeSize);
            
            [ symmetricBasisDown(indCycle,:) , phasesDownY ] = cubicSymmetries( symAppliedStates , ry, latticeSize );
            
            symmetryPhases(indCycle)  = pd * (ph_p .* phasesUpX .* ...
                phasesUpY .* ph_p_down .*phasesDownX.*phasesDownY);
            
            
        end
        
        if nSymmetries > 2
            %__________________________________________________DX
            indCycle = indCycle + nSites;
            
            [ symAppliedStates , ph_px ] = cubicSymmetries( basisStateUp , px, latticeSize);
            
            [ symAppliedStates , ph_pd ] = cubicSymmetries( symAppliedStates , pd, latticeSize );
            
            [ symAppliedStates , phasesUpX ] = cubicSymmetries( symAppliedStates , rx , latticeSize);
            
            [ symmetricBasisUp(indCycle,:) , phasesUpY ] = cubicSymmetries( symAppliedStates , ry , latticeSize);
            
            [ symAppliedStates , ph_p_down ] = cubicSymmetries( basisStateDown , px, latticeSize);
            
            [ symAppliedStates , ph_py_down ] = cubicSymmetries( symAppliedStates , pd, latticeSize);
            
            [ symAppliedStates , phasesDownX ] = cubicSymmetries( symAppliedStates , rx , latticeSize);
            
            [ symmetricBasisDown(indCycle,:) , phasesDownY ] = cubicSymmetries( symAppliedStates , ry, latticeSize );
            
            symmetryPhases(indCycle)  = pd*px * (ph_px .* ph_pd .* phasesUpX .* ...
                phasesUpY .* ph_p_down .* ph_py_down .* phasesDownX.*phasesDownY);
            
            %__________________________________________________DY
            
            indCycle = indCycle + nSites;
            
            [ symAppliedStates , ph_py ] = cubicSymmetries( basisStateUp , py, latticeSize);
            
            [ symAppliedStates , ph_pd ] = cubicSymmetries( symAppliedStates , pd, latticeSize);
            
            [ symAppliedStates , phasesUpX ] = cubicSymmetries( symAppliedStates , rx , latticeSize);
            
            [ symmetricBasisUp(indCycle,:) , phasesUpY ] = cubicSymmetries( symAppliedStates , ry , latticeSize);
            
            [ symAppliedStates , ph_p_down ] = cubicSymmetries( basisStateDown , py, latticeSize);
            
            [ symAppliedStates , ph_py_down ] = cubicSymmetries( symAppliedStates , pd, latticeSize);
            
            [ symAppliedStates , phasesDownX ] = cubicSymmetries( symAppliedStates , rx , latticeSize);
            
            [ symmetricBasisDown(indCycle,:) , phasesDownY ] = cubicSymmetries( symAppliedStates , ry, latticeSize );
            
            symmetryPhases(indCycle)  = pd*py * (ph_pd .* ph_py .* phasesUpX .* ...
                phasesUpY .* ph_p_down .* ph_py_down .* phasesDownX.*phasesDownY);
            
            
            %___________________________________________________DYX
            indCycle = indCycle + nSites;
            
            [ symAppliedStates , ph_px ] = cubicSymmetries( basisStateUp , px, latticeSize );
            
            [ symAppliedStates , ph_py ] = cubicSymmetries( symAppliedStates , py, latticeSize);
            
            [ symAppliedStates , ph_pd ] = cubicSymmetries( symAppliedStates , pd, latticeSize);
            
            [ symAppliedStates , phasesUpX ] = cubicSymmetries( symAppliedStates , rx , latticeSize);
            
            [ symmetricBasisUp(indCycle,:) , phasesUpY ] = cubicSymmetries( symAppliedStates , ry, latticeSize );
            
            [ symAppliedStates , ph_p_down ] = cubicSymmetries( basisStateDown , px, latticeSize);
            
            [ symAppliedStates , ph_py_down ] = cubicSymmetries( symAppliedStates , py, latticeSize);
            
            [ symAppliedStates , ph_pd_down ] = cubicSymmetries( symAppliedStates , pd, latticeSize);
            
            [ symAppliedStates , phasesDownX ] = cubicSymmetries( symAppliedStates , rx, latticeSize );
            
            [ symmetricBasisDown(indCycle,:) , phasesDownY ] = cubicSymmetries( symAppliedStates , ry , latticeSize);
            
            symmetryPhases(indCycle)  = pd*py*px * (ph_px .* ph_py .* ph_pd .* phasesUpX .* ...
                phasesUpY .* ph_p_down .* ph_py_down .* ph_pd_down .* phasesDownX.*phasesDownY);
            
        end
        %___________________________________________________E
        if nameFlag_pe
            indCycle = indCycle + nSites;
            
            [ symAppliedStates , ph_e ] = cubicSymmetries( basisStateUp , pe, latticeSize);
            
            [ symAppliedStates , phasesUpX ] = cubicSymmetries( symAppliedStates , rx, latticeSize);
            
            [ symmetricBasisUp(indCycle,:) , phasesUpY ] = cubicSymmetries( symAppliedStates , ry, latticeSize );
            
            [ symAppliedStates , ph_e_down ] = cubicSymmetries( basisStateDown , pe, latticeSize);
            
            [ symAppliedStates , phasesDownX ] = cubicSymmetries( symAppliedStates , rx , latticeSize);
            
            [ symmetricBasisDown(indCycle,:) , phasesDownY ] = cubicSymmetries( symAppliedStates , ry, latticeSize );
            
            symmetryPhases(indCycle)  = pe * (ph_e .* phasesUpX .* phasesUpY ...
                .* ph_e_down .*phasesDownX.*phasesDownY);
            
        end
    end
end

% combine fermionic transl and symmetryphases/eigenvalue
phasePrefactors = symmetryPhases .* exponentialPhases;

% get indices of occupied sites for up and down spin and bring in form
indSites = repmat(indSites',1,sizeSymmetryGroup);

indOccupiedUp = indSites(logical(symmetricBasisUp'));
indOccupiedDown = indSites(logical(symmetricBasisDown'));

indOccupiedUp = reshape(indOccupiedUp, nOccupiedUp, sizeSymmetryGroup)';
indOccupiedDown = reshape(indOccupiedDown, nOccupiedDown, sizeSymmetryGroup)';


