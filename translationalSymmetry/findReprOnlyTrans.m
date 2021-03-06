function [basisRepr, symOpInvariants, index2Repr, symOp2Repr] = ...
    findReprOnlyTrans(basisStates, latticeSize)
% findReprOnlyTrans calculates the representatives for translational
%  symmetry for one and two dimensions
%
% Input:        basisStates      basis in binary Form
%               nDims            lattice parameters
%
% Output:       basisRepr        representatives of trans. cycles
%               symOpInvariants  matrix containing symmetry operations 
%                                leaving each UP spin representatives 
%                                invariant. contains symOps like: 
%                                kx ky px py pd pe as eigenvalues and
%                                fermionic and exp. phase in 7th, 8th col.
%               index2Repr       linking between states and representatives
%                                sign of integer is the corresponding phase
%                                factor!
%               symOp2Repr       symmetry operations to get to repr.
%
% NOTE: 19.05.13: possibility to save alot of calculation time:
%   maybe really save "matrix" <i|T^j|k> for all basis sets in function
%   findRepresentants1D.m. than one wouldnt have to do anything in function
%   combineBasisStates1D, just link the respective basises and matrices.
%   Tj matrix also usable in hamiltonian calculation later on!
% NOTE: 25.06.13: hange secon input to lattice parameters to also include 
% rectangular lattices
% NOTE: 10.07.13: change latticeSize = [nSitesY;nSitesX] like matlab row
% column style
% 
%------------------------SVN Info------------------------------------------
% $Rev:: 63                                     $: Revision of last commit 
% $Author:: dobrautz                            $: Author of last commit   
% $Date:: 2013-06-13 09:05:27 +0200 (Don, 13. J#$: Date of last commit     
% -------------------------------------------------------------------------

%% container variables
% size of Basis and number of Sites
[nBasisStates,nSites] = size(basisStates);
% Dimension of lattice
[latticeDim,~] = size(latticeSize); 
% 1 dimensional lattice sizes: change to rect. 25.06.13: 
nSitesY = latticeSize(1);
% max.size of available states (number of sites = number of trans. syms!)
nSymmetricStates = nBasisStates*nSites;
% index for translation power identification
indTransPower = 1:nSites;
% container variable for new states(symmetry applied to them)
symmetricBasis = zeros(nSymmetricStates,nSites);
% cont. var. for respective phase factor due to symmetry exchanges(fermions)
symmetryPhases = ones(nSymmetricStates,1);
% umrechnung bin2dez
bin2dez = nSites-1:-1:0;
bin2dez = (2.^bin2dez)';
%--------------------------------------------------------------------------
% integer values of repr. for checking if we already had integer lowest repr. in cycle
integerRepr = zeros(nBasisStates,1);
% index for linking basis B to repr. of cycles, sign of this is phase fact.
index2Repr = zeros(nBasisStates,1);
% corresponding symmetry operation to get to cycle repr. (all zeros if non.comp)
symOp2Repr = zeros(nBasisStates,1);
% binary form of repr.
basisRepr = zeros(nBasisStates,nSites);
%number of found comp.repr.
nRepr = 0;
% container for phases and translation powers of invariant representatives
% the columns position of non zero elements determines the combination of
% translational powers rx, and ry: columnPos = rx*Ly + ry + 1
symOpInvariants = zeros(nBasisStates,nSites);

%% all translations on basis for different dimensions
if latticeDim == 1
    for rx1D = 0:nSites-1
        
        % application of Translational symmetry T :
        % index transformation for diff R values
        ind = rx1D + 1;
        indCycle = ind:(nSites):(nSymmetricStates-nSites+ind);
        % application on TransSymx
        [ symmetricBasis(indCycle,:) , symmetryPhases(indCycle) ] = ...
            cubicSymmetries( basisStates, rx1D, latticeSize);
        
    end
    
elseif latticeDim == 2
    nSitesX = latticeSize(2);
    for rx = 0:nSitesX-1
        for ry = 0:nSitesY-1
            
            % application of Translational symmetry T :
            % index transformation for diff R values
            ind = nSitesY*rx + ry + 1;
            % index within Symmetry
            indCycle = ind:(nSites):(nSymmetricStates-nSites+ind);
            % application on Translational Symmetry X
            [symAppliedStates, phasesX ] = cubicSymmetries(basisStates, rx, latticeSize);
            % Ty
            [symmetricBasis(indCycle,:), phasesY ] = ...
                cubicSymmetries( symAppliedStates , ry, latticeSize);
            % save fermionic phases
            symmetryPhases(indCycle) = phasesX .* phasesY;
            
        end
    end
end

%% loop over basis set
for iBasis = 0:nBasisStates-1
    
    % index to pick out one spec. symmetry operation
    indSymmetry = iBasis*nSites+1:(iBasis+1)*nSites;
    % pick binary form of one symmetry op. and calculate integer
    specSymBasis = symmetricBasis(indSymmetry,:);
    specInteger = specSymBasis*bin2dez;
    specPhases = symmetryPhases(indSymmetry);
    % find unique integers
    [uniqueInteger,indUniqueInteger,conversionIndex2UniqueInteger] = ...
        unique(specInteger);
    % look up if representant already calculated
    alreadyRepr = find(uniqueInteger(1) == integerRepr, 1);
    
    % check and update list of repr.
    if isempty(alreadyRepr)
        % increase index for every found comp. repr.
        nRepr = nRepr+1;
        % integer value of repr. (needed in the loop)
        integerRepr(nRepr) = uniqueInteger(1);
        % binary repr. (output) if its a repr. its always first element
        % since list is sorted!:
        basisRepr(nRepr,:) = specSymBasis(1,:);
        % mask for same element as starting state
        sameElementAsFirst = conversionIndex2UniqueInteger == 1;
        % store phases and translation value for invariant states
        % column position of non zero elements determines the combination
        % of translational powers rx, and ry: columnPos = rx*Ly + ry + 1
        symOpInvariants(nRepr,indTransPower(sameElementAsFirst)) = ...
            specPhases(sameElementAsFirst);
        % save index for hash table connecting basis states and repr.
        index2Repr(iBasis+1) = nRepr;
        %-------------------------note 13.05.13---------------------------
        % position in this matrix above implicetly defines the Power of the
        % translational operator :columnPos = rx*Ly + ry + 1 for 2D
        % from the amount of elements one can also conclude the number of
        % unique elements for this cycle = nSites/sum(~=0)
        %------------------------------------------------------------------
        %-------------------------note 13.05.13---------------------------
        % position in this matrix above implicetly defines the Power of the
        % translational operator T^j = T^(columnPosition-1)
        % from the amount of elements one can also conclude the number of
        % unique elements for this cycle = nSites/sum(~=0)
        %-----------------------------------------------------------------
    else
        %position of corresponding repr. times phase factor
        index2Repr(iBasis+1) = alreadyRepr*specPhases(indUniqueInteger(1));
        % symmetry operations needed to get there, for 1D easy: position in
        % symmetry group-1, for 2D only translation too : rx*Ly + ry
        symOp2Repr(iBasis+1) = indUniqueInteger(1)-1;
        
    end
    
end
% cut not used elements of container
basisRepr(nRepr+1:end,:) = [];
symOpInvariants(nRepr+1:end,:) = [];
