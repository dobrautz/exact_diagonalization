function [ basisRepr, symOpInvariants, index2Repr, symOp2Repr ] = ...
    PointGroupFindRepr( basisStates, whichSymmetry, kValue, latticeSize, varargin )
%POINTGROUPFINDREPR two dimensional representative finding function
% including point group symmetries, also calculates symmetry operation
% leaving representatives invariant and hash list connecting whole basis
% with representatives 
% 
% Input:    basisStates         binary basis state represetation
%           whichSymmetry       string flag of symmetry operation
%           kValue              k-vector values [kx;ky]
%           varargin            variable sized input of symmetry
%                               eigenvalues
% Output:   basisRepr           binary representation of found
%                               representatives
%           symOpInvariants     symmetry operations leaving representatives
%                               invariant: cell array!
%                               cols 1:6 : symmetries with eigenvalues
%                               col 7:     fermionic phase
%                               col 8:     expontial phase
%           index2Repr          hash list connecting all UP spin basis
%                               states to corresponding representatives and
%                               also includes fermionic phase as the SIGN
%           symOp2Repr          corresponding symmetry operation bringing
%                               the state to its representative
% 
% NOTE: 24.05.13: for PG inclusion old elegant way of saving symmetry 
% operations doesnt work anymore, since there are too many possible 
% symmetry combinations -> do CELL structure again
% NOTE: 24.05.13: first draft: copied together from previous versions of
% this functions; most similar to old function applyCubicLatticeSymm 
% NOTE: 12.06.13: almost final form of this function:
% * with LaFlorience combination method i dont need DOWN spin representatives
% anymore: so changed output callsof this function. 
% * for faster norm calculation already in combination2Hubbard function i
% also save the exponential phases in the symOpInvariants matrix!
% * changed to pre calculation of symmetryOperation matrix and exponential
% phases( possible error-prone)
% 
%------------------------SVN Info------------------------------------------
% $Rev:: 63                                     $: Revision of last commit
% $Author:: dobrautz                            $: Author of last commit
% $Date:: 2013-06-13 09:05:27 +0200 (Don, 13. J#$: Date of last commit
% -------------------------------------------------------------------------

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
% size of Basis and number of Sites (1D and 2D)
[nBasisStates, nSites] = size(basisStates);
nSitesY = latticeSize(1);
nSitesX = latticeSize(2);
% size of symmetry group
sizeSymmetryGroup = combSymmetries*nSites;    
% max.size of available states
nSymmetricStates = sizeSymmetryGroup * nBasisStates;
% container variable for new states(symmetry applioed to them)
symmetricBasis = zeros(nSymmetricStates,nSites);
% cont. var. for respective phase factor due to symmetry exchanges(fermions)
symmetryPhases = ones(nSymmetricStates,1);
% umrechnung bin2dez
bin2dez = nSites-1:-1:0;
bin2dez = (2.^bin2dez)';
% -------------------------------------------------------------------------
% integer values of repr. for checking if we already had integer lowest repr. in cycle
integerRepr = zeros(nBasisStates,1);
% index for linking basis B to repr. of cycles, sign of this is phase fact.
index2Repr = zeros(nBasisStates,1);
% corresponding symmetry operation to get to cycle repr. (all zeros if non.comp)
symOp2Repr = zeros(nBasisStates,6);
% binary form of repr.
basisRepr = zeros(nBasisStates,nSites);
%number of found comp.repr.
nRepr = 0;
% container for phases and translation powers of invariant representatives
% way of saving it: for each representative: 
% kx, ky, px, py, pd, pe, phase, expPhase
symOpInvariants = cell(nBasisStates,1);

% cont. var. for symmetry operations applied
% storage in R: kx , ky , px , py , pd , pe
symmetryOperations = zeros(nSymmetricStates,6);
% keyboard
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
exponentialPhases = repmat(exponentialPhases,nBasisStates,1);
%--------------------------------------------------------------------------

%% all translations and corresponding point group syms on basis
% loops over all possible R values (square lattice!)
for rx = 0:nSitesX-1
    for ry = 0:nSitesY-1
        
        % application of Translational symmetry T :
        % index transformation for diff R values
        ind = nSitesY*rx + ry + 1;
        % index within Symmetry
        indCycle = ind : sizeSymmetryGroup : (nSymmetricStates - sizeSymmetryGroup + ind);
        
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
        
        % application of Translation Operator : X first
        [ symAppliedStates , ph_x ] = cubicSymmetries( basisStates , rx, latticeSize ,'rx'); 
        % then Y translation
        [ symmetricBasis(indCycle,:) , ph_y ] = cubicSymmetries( symAppliedStates , ry, latticeSize ,'ry');
        % save fermionic phases
        symmetryPhases(indCycle) = ph_x .* ph_y; 
        
         %applied symmetry operation(sign convention)Remember this!
        symmetryOperations(indCycle,1) = rx; 
        symmetryOperations(indCycle,2) = ry; 
       
        % bis hier immer
        
        % now decide which point group is to be applied: 
        if nameFlag_px || nSymmetries > 1
            %application of Translational symmetry T and Px
            indCycle = indCycle + nSites; %increase ind with number of translational syms=x
            
            [ symAppliedStates , ph_p ] = cubicSymmetries(basisStates , px, latticeSize,'px');
            
            [ symAppliedStates , ph_x ] = cubicSymmetries( symAppliedStates , rx, latticeSize ,'rx');
            
            [ symmetricBasis(indCycle,:) , ph_y ] = cubicSymmetries( symAppliedStates , ry, latticeSize,'ry' );
            
            symmetryPhases(indCycle)  = ph_p .* ph_x .* ph_y;

            symmetryOperations(indCycle,1) = rx ;
            symmetryOperations(indCycle,2) = ry; 
            symmetryOperations(indCycle,3) = px;
            
        end
        
        if nameFlag_py || nSymmetries > 1
            %TPy
            indCycle = indCycle + nSites;
            
            [ symAppliedStates , ph_p ] = cubicSymmetries( basisStates , py, latticeSize, 'py');
            
            [ symAppliedStates , ph_x ] = cubicSymmetries( symAppliedStates , rx, latticeSize ,'rx');
            
            [ symmetricBasis(indCycle,:) , ph_y ] = cubicSymmetries( symAppliedStates , ry, latticeSize ,'ry');
            
            symmetryPhases(indCycle)  = ph_p .* ph_x .* ph_y;

            symmetryOperations(indCycle,1) = rx ;
            symmetryOperations(indCycle,2) = ry; 
            symmetryOperations(indCycle,4) = py;
        end
        
        if nSymmetries > 1
            %_________________________________________XY
            
            indCycle = indCycle + nSites;
            
            [ symAppliedStates , ph_px ] = cubicSymmetries( basisStates , px , latticeSize,'px');
            
            [ symAppliedStates , ph_py ] = cubicSymmetries( symAppliedStates , py, latticeSize , 'py');
            
            [ symAppliedStates , ph_x ] = cubicSymmetries( symAppliedStates , rx, latticeSize ,'rx');
            
            [ symmetricBasis(indCycle,:) , ph_y ] = cubicSymmetries( symAppliedStates , ry , latticeSize,'ry');
            
            symmetryPhases(indCycle)  = ph_px .* ph_py .* ph_x .* ph_y;
            
            symmetryOperations(indCycle,1) = rx ; 
            symmetryOperations(indCycle,2) = ry; 
            symmetryOperations(indCycle,3) = px;   
            symmetryOperations(indCycle,4) = py;
        end
        
        if nameFlag_pd || nSymmetries > 2
            %__________________________________D
            
            indCycle = indCycle + nSites;
            
            [ symAppliedStates , ph_p ] = cubicSymmetries( basisStates , pd, latticeSize ,'pd');
            
            [ symAppliedStates , ph_x ] = cubicSymmetries( symAppliedStates , rx, latticeSize ,'rx');
            
            [ symmetricBasis(indCycle,:) , ph_y ] = cubicSymmetries( symAppliedStates , ry, latticeSize ,'ry');
            
            symmetryPhases(indCycle)  = ph_p .* ph_x .* ph_y;

            symmetryOperations(indCycle,1) = rx ; 
            symmetryOperations(indCycle,2) = ry; 
            symmetryOperations(indCycle,5) = pd;
        end
        
        if nSymmetries > 2
            %__________________________________________________DX
            indCycle = indCycle + nSites;
            
            [ symAppliedStates , ph_px ] = cubicSymmetries( basisStates , px , latticeSize,'px');
            
            [ symAppliedStates , ph_pd ] = cubicSymmetries( symAppliedStates , pd , latticeSize, 'pd');
            
            [ symAppliedStates , ph_x ] = cubicSymmetries( symAppliedStates , rx , latticeSize,'rx');
            
            [ symmetricBasis(indCycle,:) , ph_y ] = cubicSymmetries( symAppliedStates , ry, latticeSize ,'ry');
            
            symmetryPhases(indCycle)  = ph_px .* ph_pd .* ph_x .* ph_y;
            
            symmetryOperations(indCycle,1) = rx ; 
            symmetryOperations(indCycle,2) = ry; 
            symmetryOperations(indCycle,3) = px;   
            symmetryOperations(indCycle,5) = pd;
            %__________________________________________________DY
            
            indCycle = indCycle + nSites;
            
            [ symAppliedStates , ph_py ] = cubicSymmetries( basisStates , py , latticeSize,'py');
            
            [ symAppliedStates , ph_pd ] = cubicSymmetries( symAppliedStates , pd , latticeSize,'pd');
            
            [ symAppliedStates , ph_x ] = cubicSymmetries( symAppliedStates , rx , latticeSize,'rx');
            
            [ symmetricBasis(indCycle,:) , ph_y ] = cubicSymmetries( symAppliedStates , ry , latticeSize,'ry' );
            
            symmetryPhases(indCycle)  = ph_py .* ph_pd .* ph_x .* ph_y;
            
            symmetryOperations(indCycle,1) = rx ; 
            symmetryOperations(indCycle,2) = ry; 
            symmetryOperations(indCycle,4) = py;   
            symmetryOperations(indCycle,5) = pd;
            %___________________________________________________DYX
            indCycle = indCycle + nSites;
            
            [ symAppliedStates , ph_px ] = cubicSymmetries( basisStates , px , latticeSize,'px');
            
            [ symAppliedStates , ph_py ] = cubicSymmetries( symAppliedStates , py, latticeSize,'py');
            
            [ symAppliedStates , ph_pd ] = cubicSymmetries( symAppliedStates , pd, latticeSize ,'pd');
            
            [ symAppliedStates , ph_x ] = cubicSymmetries( symAppliedStates , rx , latticeSize,'rx');
            
            [ symmetricBasis(indCycle,:) , ph_y ] = cubicSymmetries( symAppliedStates , ry , latticeSize,'ry');
            
            symmetryPhases(indCycle)  = ph_px .* ph_py .* ph_pd .* ph_x .* ph_y;
            
            symmetryOperations(indCycle,1) = rx ; 
            symmetryOperations(indCycle,2) = ry; 
            symmetryOperations(indCycle,3) = px;  
            symmetryOperations(indCycle,4) = py; 
            symmetryOperations(indCycle,5) = pd;
        end
        %___________________________________________________E
        if nameFlag_pe
            indCycle = indCycle + nSites;
            
            [ symAppliedStates , ph_e ] = cubicSymmetries( basisStates , pe , latticeSize,'pe');
            
            [ symAppliedStates , ph_x ] = cubicSymmetries( symAppliedStates , rx, latticeSize, 'rx' );
            
            [ symmetricBasis(indCycle,:) , ph_y ] = cubicSymmetries( symAppliedStates , ry, latticeSize, 'ry' );
            
            symmetryPhases(indCycle)  = ph_e .* ph_x .* ph_y;
            
            symmetryOperations(indCycle,1) = rx ; 
            symmetryOperations(indCycle,2) = ry; 
            symmetryOperations(indCycle,6) = pe;
        end
        
    end
end

%% loop over basis set
for iBasis = 0:nBasisStates-1      
    
    % index to pick out one spec. symmetry operation
    indSymmetry = (iBasis*sizeSymmetryGroup+1):(iBasis+1)*sizeSymmetryGroup;          
    % pick binary form of one symmetry op.
    specSymBasis = symmetricBasis(indSymmetry,:);
    % corresponding phases, symmetry operations and integer values 
    specPhases = symmetryPhases(indSymmetry);
    specExpPhases = exponentialPhases(indSymmetry);
    specSymOperation = symmetryOperations(indSymmetry,:);              
    specInteger = specSymBasis*bin2dez;      
    % find unique integers
    [uniqueInteger, indUniqueInteger, conversionIndex2UniqueInteger] = ...
        unique(specInteger,'first');
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
        % build together symOpInvariance matrix:
        specSymOpInv = zeros(sum(sameElementAsFirst),8);
        specSymOpInv(:,1:6) = specSymOperation(sameElementAsFirst,:);
        specSymOpInv(:,7) = specPhases(sameElementAsFirst);
        specSymOpInv(:,8) = specExpPhases(sameElementAsFirst);
        % save it in cell:
        symOpInvariants{nRepr} = specSymOpInv;
        % save index for hash table connecting basis states and repr.
        index2Repr(iBasis+1) = nRepr;
        
    else
        % position of corresponding repr. times phase factor
        index2Repr(iBasis+1) = alreadyRepr*specPhases(indUniqueInteger(1));
        % symmetry operations needed to get there, for 1D easy: position in
        % symmetry group-1, for 2D only translation too : rx*Ly + ry
        % for 2D with PG save it now the old way: kx, ky, px, py, pd, pe
        symOp2Repr(iBasis+1,:) = specSymOperation(indUniqueInteger(1),:);

    end
    
end
% cut not used elements of container
basisRepr(nRepr+1:end,:) = [];
symOpInvariants(nRepr+1:end,:) = [];


% NOTE: 12.06.13: below not needed anymore: with LaFlorience method i dont 
% need down spin representatives!
% % depending of number of output arguments(if its called for UP or DOWN spin
% % different output
% if nargout == 2 % DOWN spin case
%     
%     % ----------NOTE: 20.05.13: for 1Dim case this was sufficient:---------
%     % number of unique Elements of a cycle = factor \nu_n from fehske
%     % for application of T^j in hubbard basis combination : 
%     % |r> = |nUp> T^j |mDown> ; j = 0...min(\nu_n,\nu_m)-1 
%     % ---------------------------------------------------------------------
%     % but now the combination of X and Y translations is needed, so i also
%     % output symOpInvariants for Down Spin case
%     % column position of non zero elements determines the combination 
%     % of translational powers rx, and ry: columnPos = rx*Ly + ry + 1:
%     %----------------------------------------------------------------------
%     varargout(1) = {symOpInvariants};
%     
% else            % UP spin case
%     
%     varargout(1) = {symOpInvariants};
%     varargout(2) = {index2Repr};
%     varargout(3) = {symOp2Repr};
% end
