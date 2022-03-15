function [ compDownStatesPerRepr,compInd2ReprDown, normHubbardStates ] = ...
    PointGroupHubbardNorm( basisReprUp, compDownStatesPerRepr,compInd2ReprDown, kValue, whichSymmetry, varargin )
%POINTGROUPHUBBARDNORM sorts out non compatible representatives and calculates
% the norm of the compatible hubbard basis states for 2 dimensional case
% including point group symmetries
%
% NOTE: 24.05.13: first draft: copied together from previous versions of
% this functions; combination of twoDimCalcHubbardNorm and
% applyCubicLatticeSyms and calcRepr.
% NOTE: 12.06.13: not needed anymore: this functionality included in basis
% combination functions!
%------------------------SVN Info------------------------------------------
% $Rev:: 62                                     $: Revision of last commit
% $Author:: dobrautz                            $: Author of last commit
% $Date:: 2013-06-12 15:46:11 +0200 (Mit, 12. J#$: Date of last commit
% -------------------------------------------------------------------------


%% identify which symmetry operation is to be applied
nameFlag_px = strcmp(whichSymmetry,'px');
nameFlag_py = strcmp(whichSymmetry,'py');
nameFlag_pd = strcmp(whichSymmetry,'pd');
nameFlag_pe = strcmp(whichSymmetry,'pe');

% number of variable arguments in = number of symmetry operations!
nSymmetries = numel(varargin);

if  nSymmetries == 1 % case: one lattice sym
    
    if nameFlag_px
        px = varargin{1};
    elseif nameFlag_py
        py = varargin{1};
    elseif nameFlag_pd
        pd = varargin{1};
    elseif nameFlag_pe
        pe = varargin{1};
    end
    
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
% size of up spin Basis repr. and number of sites 2D and 1D
[nReprStatesUp,nSites] = size(basisReprUp);
nSites1D = sqrt(nSites);
% size of symmetry group
sizeSymmetryGroup = combSymmetries*nSites;
% container for norm of compatible hubbard basis states
normHubbardStates = cell(nReprStatesUp,1);
% umrechnung bin2dez
bin2dez = nSites-1:-1:0;
bin2dez = (2.^bin2dez)';

% loop over UP spin representatives
for iReprUp = 1:nReprStatesUp
    
    % pick out DOWN spin basis and calc number of respective states
    basisReprDown = compDownStatesPerRepr{iReprUp};
    [nReprStatesDown,referenceTest] = size(basisReprDown);
    % determine if specific down state representatice is reference and
    % assign corresponding reference value
    if referenceTest == 1 
        
        % same results for all non invariant sates: just reference to it
        % value in compDownStatesPerRepr is index of reference state!
        compDownStatesPerRepr{iReprUp} = basisReprDown;
        compInd2ReprDown{iReprUp} = basisReprDown;
        normHubbardStates{iReprUp} = normHubbardStates{basisReprDown};
   
    else
    % replicate up spin repr. so often
    basisReprUpRepl = repmat(basisReprUp(iReprUp,:),nReprStatesDown,1);
    % max.size of available states
    nSymmetricStates = sizeSymmetryGroup*nReprStatesDown;
    % container variable for new states(symmetry applioed to them)
    symmetricBasisUp = zeros(nSymmetricStates,nSites);
    % symmetric basis Down:
    symmetricBasisDown = zeros(nSymmetricStates,nSites);
    % cont. var. for respective phase factor due to symmetry exchanges(fermions)
    symmetryPhases = ones(nSymmetricStates,1);
    % exp. phases pre-determined from lattice size and symmetry size:
    rVectorX = 0:nSites1D-1;
    rVectorY = 0:nSites1D-1;
    [rVectorY, rVectorX] = ndgrid(rVectorX, rVectorY);
    rVectorX = rVectorX(:);
    rVectorY = rVectorY(:);
    exponentialPhases = exp(1i*(kValue(1)*rVectorX + kValue(2)*rVectorY));
    % surpress real and imaginary parts of this below treshhold
    exponentialPhases = real(exponentialPhases).*(abs(real(exponentialPhases))>10^-5) + ...
        1i*imag(exponentialPhases).*(abs(imag(exponentialPhases))>10^-5);
    % repmat it to whole symmetry basis ( actually not needed but for
    % uniform indication useful
    exponentialPhases = repmat(exponentialPhases,combSymmetries,1);
    exponentialPhases = repmat(exponentialPhases,nReprStatesDown,1);
    
    %% all translations and corresponding point group syms on basis
    % loops over all possible R values (square lattice!)
    for rx = 0:nSites1D-1
        for ry = 0:nSites1D-1
            
            % application of Translational symmetry T :
            % index transformation for diff R values
            ind = nSites1D*rx + ry + 1;
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
            
            
            % application of symmetry operation on UP spin part: X first
            [ symAppliedStates , phasesUpX ] = cubicSymmetries( basisReprUpRepl , rx );
            % then Y translation
            [ symmetricBasisUp(indCycle,:) , phasesUpY ] = cubicSymmetries( symAppliedStates , ry );
            % same with DOWN SPIN part:
            [ symAppliedStates , phasesDownX ] = cubicSymmetries( basisReprDown , rx );
            % then Y translation
            [ symmetricBasisDown(indCycle,:) , phasesDownY ] = cubicSymmetries( symAppliedStates , ry );
            
            % save fermionic phases
            symmetryPhases(indCycle) = (phasesUpX .* phasesUpY .* phasesDownX .* phasesDownY) ;
 
            % bis hier immer
            
            % now decide which point group is to be applied:
            if nameFlag_px || nSymmetries > 1
                %application of Translational symmetry T and Px
                indCycle = indCycle + nSites; %increase ind with number of translational syms=x
                
                [ symAppliedStates , ph_p ] = cubicSymmetries( basisReprUpRepl , px);
                
                [ symAppliedStates , phasesUpX ] = cubicSymmetries( symAppliedStates , rx );
                
                [ symmetricBasisUp(indCycle,:) , phasesUpY ] = cubicSymmetries( symAppliedStates , ry );
                
                [ symAppliedStates , ph_p_down ] = cubicSymmetries( basisReprDown , px);
                
                [ symAppliedStates , phasesDownX ] = cubicSymmetries( symAppliedStates , rx );
                
                [ symmetricBasisDown(indCycle,:) , phasesDownY ] = cubicSymmetries( symAppliedStates , ry );
                
                symmetryPhases(indCycle)  = px * (ph_p .* phasesUpX .* phasesUpY .* ...
                    ph_p_down .*phasesDownX.*phasesDownY);
                
                
            end
            
            if nameFlag_py || nSymmetries > 1
                %TPy
                indCycle = indCycle + nSites;
                 
                [ symAppliedStates , ph_p ] = cubicSymmetries( basisReprUpRepl , py);
                
                [ symAppliedStates , phasesUpX ] = cubicSymmetries( symAppliedStates , rx );
                
                [ symmetricBasisUp(indCycle,:) , phasesUpY ] = cubicSymmetries( symAppliedStates , ry );
                
                [ symAppliedStates , ph_p_down ] = cubicSymmetries( basisReprDown , py);
                
                [ symAppliedStates , phasesDownX ] = cubicSymmetries( symAppliedStates , rx );
                
                [ symmetricBasisDown(indCycle,:) , phasesDownY ] = cubicSymmetries( symAppliedStates , ry );
                
                symmetryPhases(indCycle)  = py * (ph_p .* phasesUpX .* phasesUpY...
                    .* ph_p_down .*phasesDownX.*phasesDownY);
                
                
            end
            
            if nSymmetries > 1
                %_________________________________________XY
                
                indCycle = indCycle + nSites;
                
                [ symAppliedStates , ph_px ] = cubicSymmetries( basisReprUpRepl , px);
                
                [ symAppliedStates , ph_py ] = cubicSymmetries( symAppliedStates , py);
                
                [ symAppliedStates , phasesUpX ] = cubicSymmetries( symAppliedStates , rx );
                
                [ symmetricBasisUp(indCycle,:) , phasesUpY ] = cubicSymmetries( symAppliedStates , ry );
                
                [ symAppliedStates , ph_p_down ] = cubicSymmetries( basisReprDown , px);
                
                [ symAppliedStates , ph_py_down ] = cubicSymmetries( symAppliedStates , py);
                
                [ symAppliedStates , phasesDownX ] = cubicSymmetries( symAppliedStates , rx );
                
                [ symmetricBasisDown(indCycle,:) , phasesDownY ] = cubicSymmetries( symAppliedStates , ry );
                
                symmetryPhases(indCycle)  = px*py * (ph_px .* ph_py .* phasesUpX .* ...
                    phasesUpY .* ph_p_down .* ph_py_down .* phasesDownX.*phasesDownY);
                
                
            end
            
            if nameFlag_pd || nSymmetries > 2
                %__________________________________D
                
                indCycle = indCycle + nSites;
                
                [ symAppliedStates , ph_p ] = cubicSymmetries( basisReprUpRepl , pd );
                
                [ symAppliedStates , phasesUpX ] = cubicSymmetries( symAppliedStates , rx );
                
                [ symmetricBasisUp(indCycle,:) , phasesUpY ] = cubicSymmetries( symAppliedStates , ry );
                
                [ symAppliedStates , ph_p_down ] = cubicSymmetries( basisReprDown , pd);
                
                [ symAppliedStates , phasesDownX ] = cubicSymmetries( symAppliedStates , rx );
                
                [ symmetricBasisDown(indCycle,:) , phasesDownY ] = cubicSymmetries( symAppliedStates , ry );
                
                symmetryPhases(indCycle)  = pd * (ph_p .* phasesUpX .* ...
                    phasesUpY .* ph_p_down .*phasesDownX.*phasesDownY);

                
            end
            
            if nSymmetries > 2
                %__________________________________________________DX
                indCycle = indCycle + nSites;

                [ symAppliedStates , ph_px ] = cubicSymmetries( basisReprUpRepl , px);
                
                [ symAppliedStates , ph_pd ] = cubicSymmetries( symAppliedStates , pd );
                
                [ symAppliedStates , phasesUpX ] = cubicSymmetries( symAppliedStates , rx );
                
                [ symmetricBasisUp(indCycle,:) , phasesUpY ] = cubicSymmetries( symAppliedStates , ry );
                
                [ symAppliedStates , ph_p_down ] = cubicSymmetries( basisReprDown , px);
                
                [ symAppliedStates , ph_py_down ] = cubicSymmetries( symAppliedStates , pd);
                
                [ symAppliedStates , phasesDownX ] = cubicSymmetries( symAppliedStates , rx );
                
                [ symmetricBasisDown(indCycle,:) , phasesDownY ] = cubicSymmetries( symAppliedStates , ry );
                
                symmetryPhases(indCycle)  = pd*px * (ph_px .* ph_pd .* phasesUpX .* ...
                    phasesUpY .* ph_p_down .* ph_py_down .* phasesDownX.*phasesDownY);
                
                %__________________________________________________DY
                
                indCycle = indCycle + nSites;
                
                [ symAppliedStates , ph_py ] = cubicSymmetries( basisReprUpRepl , py);
                
                [ symAppliedStates , ph_pd ] = cubicSymmetries( symAppliedStates , pd);
                
                [ symAppliedStates , phasesUpX ] = cubicSymmetries( symAppliedStates , rx );
                
                [ symmetricBasisUp(indCycle,:) , phasesUpY ] = cubicSymmetries( symAppliedStates , ry );
                
                [ symAppliedStates , ph_p_down ] = cubicSymmetries( basisReprDown , py);
                
                [ symAppliedStates , ph_py_down ] = cubicSymmetries( symAppliedStates , pd);
                
                [ symAppliedStates , phasesDownX ] = cubicSymmetries( symAppliedStates , rx );
                
                [ symmetricBasisDown(indCycle,:) , phasesDownY ] = cubicSymmetries( symAppliedStates , ry );
                
                symmetryPhases(indCycle)  = pd*py * (ph_pd .* ph_py .* phasesUpX .* ...
                    phasesUpY .* ph_p_down .* ph_py_down .* phasesDownX.*phasesDownY);
                
                
                %___________________________________________________DYX
                indCycle = indCycle + nSites;
                
                [ symAppliedStates , ph_px ] = cubicSymmetries( basisReprUpRepl , px );
                
                [ symAppliedStates , ph_py ] = cubicSymmetries( symAppliedStates , py);
                
                [ symAppliedStates , ph_pd ] = cubicSymmetries( symAppliedStates , pd);
                
                [ symAppliedStates , phasesUpX ] = cubicSymmetries( symAppliedStates , rx );
                
                [ symmetricBasisUp(indCycle,:) , phasesUpY ] = cubicSymmetries( symAppliedStates , ry );
                
                [ symAppliedStates , ph_p_down ] = cubicSymmetries( basisReprDown , px);
                
                [ symAppliedStates , ph_py_down ] = cubicSymmetries( symAppliedStates , py);
                
                [ symAppliedStates , ph_pd_down ] = cubicSymmetries( symAppliedStates , pd);
                
                [ symAppliedStates , phasesDownX ] = cubicSymmetries( symAppliedStates , rx );
                
                [ symmetricBasisDown(indCycle,:) , phasesDownY ] = cubicSymmetries( symAppliedStates , ry );
                
                symmetryPhases(indCycle)  = pd*py*px * (ph_px .* ph_py .* ph_pd .* phasesUpX .* ...
                    phasesUpY .* ph_p_down .* ph_py_down .* ph_pd_down .* phasesDownX.*phasesDownY);
                
            end
            %___________________________________________________E
            if nameFlag_pe
                indCycle = indCycle + nSites;
               
                [ symAppliedStates , ph_e ] = cubicSymmetries( basisReprUpRepl , pe);
                
                [ symAppliedStates , phasesUpX ] = cubicSymmetries( symAppliedStates , rx );
                
                [ symmetricBasisUp(indCycle,:) , phasesUpY ] = cubicSymmetries( symAppliedStates , ry );
                
                [ symAppliedStates , ph_e_down ] = cubicSymmetries( basisReprDown , pe);
                
                [ symAppliedStates , phasesDownX ] = cubicSymmetries( symAppliedStates , rx );
                
                [ symmetricBasisDown(indCycle,:) , phasesDownY ] = cubicSymmetries( symAppliedStates , ry );
                
                symmetryPhases(indCycle)  = pe * (ph_e .* phasesUpX .* phasesUpY ...
                    .* ph_e_down .*phasesDownX.*phasesDownY);

            end
        end
    end
    
    %% loop over basis set
    indCompState = 1; % counting index for found comp. sates
    for iBasis = 0:nReprStatesDown-1
        
        % index to pick out one spec. symmetry operation
        indSymmetry = ( iBasis*sizeSymmetryGroup + 1):( (iBasis+1) * sizeSymmetryGroup); 
        
        % pick binary form of one symmetry op.
        specSymBasisUp = symmetricBasisUp(indSymmetry,:);
        specSymBasisDown = symmetricBasisDown(indSymmetry,:);
        % corresponding phases, symmetry operations and integer values
        specPhases = symmetryPhases(indSymmetry) .* exponentialPhases(indSymmetry);
        specIntegerUp = specSymBasisUp*bin2dez;
        specIntegerDown = specSymBasisDown*bin2dez;
        % combine integer for MSB UP spin and LSB DOWN spin part
        specInteger = (2^nSites)*specIntegerUp + specIntegerDown;
        % here i only need the same elements as first one: no need of
        % unique fct. here and in 1Dim version of this
        % calulate unique elements and respective index
        sameElementAsFirst = specInteger == specInteger(1);
        
        % ---------------- note: 13.05.2013--------------------------------
        % in SANDVIK paper ED techniques:
        % Norm = (number of unique elements) * |sum(phases(sameAsFirst)|
        %
        % [uniqueInteger,~,conversionIndex2UniqueInteger] = ...
        %       unique(specInteger,'stable');
        % mask for same element as starting state |UP>|DOWN>
        % sameElementAsFirst = conversionIndex2UniqueInteger == 1;
        % % # of unique elements of unique cycle
        % nUniqueElements = numel(uniqueInteger); 
        % sumPhases = sumPhases*conj(sumPhases);
        % combine it for norm:
        % normSpecHubbardState = nUniqueElements*sumPhases;
        % ----------------- but in FEHSKE PAPER----------------------------
        % only sum(phases(sameAsFirst))/nSites
        % DO IT AS IN FEHSKE for now
        %--------------------end note 13.05.2013---------------------------
        
        % sum of phases of same elements
        sumPhases = sum(specPhases(sameElementAsFirst));
        % i think i can take abs value
        normSpecHubbardState = abs(sumPhases)/nSites; 
        
        % if state is compatible(norm~=0) save norm in container
        if normSpecHubbardState > 10^-10
            normHubbardStates{iReprUp}(indCompState,1) = normSpecHubbardState;
            indCompState = indCompState + 1;
            
        else % if state is incompatible (norm == 0) delete corr. repr.
            compDownStatesPerRepr{iReprUp}(indCompState,:) = [];
            compInd2ReprDown{iReprUp}(indCompState) = [];
        end
    end
    end
end


    
