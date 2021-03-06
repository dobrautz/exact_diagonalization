function [ compDownStatesPerRepr, compInd2ReprDown, normHubbardStates ] = ...
    calcHubbardNormOnlyTrans(basisReprUp, compDownStatesPerRepr, compInd2ReprDown, kValue)
%CALCHUBBARDNORMONLYTRANS sorts out non compatible representatives and 
% calculates the norm of the compatible hubbard basis states for one and 
% two dimensions and only translational symmetry
% 
% Input:    basisReprUp             binary representation of up spin
%                                   representatives
%           compDownStatesPerRepr   whole down state representatives per up
%                                   spin representatives. this function
%                                   sorts out the imcompatible ones of this
%                                   list
%           compInd2ReprDown        index linking to whole down spin basis
%           kValue                  k value for which comp. repr. are
%                                   calulated
% 
% Output:   compDownStatesPerRepr   updated compatible basis states
%           compInd2ReprDown        updated list linking to whole basis
%           normHubbard             corresponding norm of states
% NOTE: 27.05.13: added additional I/O: a linking list between the down
% spin part of the representatives and the whole down spin basis. because
% in calculation of down spin hopping part of hamiltionian i can once apply
% it to whole down spin basis and then link it to respective down spin part
% tp save time
% NOTE: 28.05.13: adapt this function to changed way of saving
% index2ReprDown and downStatesPerRepr( for duplicate values only reference
% to first appearence) 
% NOTE: 28.05.13: IDEA if the down spin part of a representative is the
% whole down spin basis -> this means there are no symmetry operations
% leaving the up spin part invariant -> there are only unique values
% created by application of the symmetry operations to up spin part 
% -> this means all values, even combined, are unique! so i would only need
% first value of phases etc of a representative to calc norm which is only 
% 1/nSites! for each down spin member!
% did this now in same way as in TwoDimCombine algortihm with references to
% first appearence!
% 
%------------------------SVN Info------------------------------------------
% $Rev::                                        $: Revision of last commit 
% $Author:: dobrautz                            $: Author of last commit   
% $Date:: 2013-06-22 19:51:27 +0200 (Sam, 22. J#$: Date of last commit     
% -------------------------------------------------------------------------

% size of up spin Basis repr. and number of sites 2D and 1D
[nReprStatesUp,nSites] = size(basisReprUp);
nSites1D = sqrt(nSites);
% container for norm of compatible hubbard basis states
normHubbardStates = cell(nReprStatesUp,1);
% umrechnung bin2dez
bin2dez = nSites-1:-1:0;
bin2dez = (2.^bin2dez)';
% number of dimensions
nDims = numel(kValue);

%% loop over UP spin representatives
for iReprUp = 1:nReprStatesUp
    
    % pick out DOWN spin basis and calc number of respective states
    basisReprDown = compDownStatesPerRepr{iReprUp};
    [nReprStatesDown,~] = size(basisReprDown);    
    % replicate up spin repr. so often
    basisReprUpRepl = repmat(basisReprUp(iReprUp,:),nReprStatesDown,1);
    % max.size of available states (number of sites = number of trans. syms!)
    nSymmetricStates = nReprStatesDown*nSites;
    % container variable for new states(symmetry applioed to them)
    symmetricBasisUp = zeros(nSymmetricStates,nSites);
    % symmetric basis DOwn:
    symmetricBasisDown = zeros(nSymmetricStates,nSites);
    % cont. var. for respective phase factor due to symmetry exchanges(fermions)
    symmetryPhases = ones(nSymmetricStates,1);
    % cont. variable for exponential phases due to translational symmetry
    exponentialPhases = zeros(nSymmetricStates,1);

    %% all translations on basis
    if nDims == 1
        
        for rx1D = 0:nSites-1
            % application of Translational symmetry T :
            % index transformation for diff R values
            ind = rx1D + 1;
            indCycle = ind:(nSites):(nSymmetricStates-nSites+ind);
            % calc. of exp phase factor for !careful for sign in exponential!
            exponentialPhases(indCycle) = exp(1i*kValue*rx1D);
            % application on TransSymx to UP and DOWN spin part:
            [symmetricBasisUp(indCycle,:), symmetryPhasesUp ] = ...
                cubicSymmetries(basisReprUpRepl, rx1D );
            [symmetricBasisDown(indCycle,:), symmetryPhasesDown ] = ...
                cubicSymmetries(basisReprDown, rx1D );
            % combine fermionic phases
            symmetryPhases(indCycle) = symmetryPhasesUp.*symmetryPhasesDown;
            
        end
        
    elseif nDims == 2
        
        for rx = 0:nSites1D-1
            for ry = 0:nSites1D-1
                % application of Translational symmetry T :
                % index transformation for diff R values
                ind = nSites1D*rx + ry + 1;
                % index within Symmetry
                indCycle = ind:(nSites):(nSymmetricStates-nSites+ind);
                % calc. of exp phase factor for !careful for sign in exponential!
                exponentialPhases(indCycle) = exp(1i*(kValue(1)*rx + kValue(2)*ry));
                % UP SPIN PART: application on Translational Symmetry X
                [ symAppliedStates , phasesUpX ] = cubicSymmetries( basisReprUpRepl , rx );
                % Ty
                [ symmetricBasisUp(indCycle,:) , phasesUpY ] = cubicSymmetries( symAppliedStates , ry );
                % DOWN SPIN PART: application on Translational Symmetry X
                [ symAppliedStates , phasesDownX ] = cubicSymmetries( basisReprDown , rx );
                % Ty
                [ symmetricBasisDown(indCycle,:) , phasesDownY ] = cubicSymmetries( symAppliedStates , ry );
                % save fermionic phases
                symmetryPhases(indCycle) = phasesUpX .* phasesUpY .* phasesDownX .* phasesDownY;
                
            end
            
        end
    end
    
    %% loop over basis set
    indCompState = 1; % counting index for found comp. sates
    for iBasis = 0:nReprStatesDown-1
        
        % index to pick out one spec. symmetry operation
        indSymmetry = iBasis*nSites+1:(iBasis+1)*nSites;
        % pick binary form of one symmetry op.
        specSymBasisUp = symmetricBasisUp(indSymmetry,:);
        specSymBasisDown = symmetricBasisDown(indSymmetry,:);
        % corresponding phases, symmetry operations and integer values
        specPhases = symmetryPhases(indSymmetry).*exponentialPhases(indSymmetry);
        specIntegerUp = specSymBasisUp*bin2dez;
        specIntegerDown = specSymBasisDown*bin2dez;
        % combine integer for MSB UP spin and LSB DOWN spin part
        specInteger = 2^nSites*specIntegerUp + specIntegerDown;
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
        if normSpecHubbardState > 10^-5
            normHubbardStates{iReprUp}(indCompState,1) = normSpecHubbardState;
            indCompState = indCompState + 1;
            
        else % if state is incompatible (norm == 0) delete corr. repr.
            compDownStatesPerRepr{iReprUp}(indCompState,:) = [];
            compInd2ReprDown{iReprUp}(indCompState) = [];
   
        end
    end    
end


