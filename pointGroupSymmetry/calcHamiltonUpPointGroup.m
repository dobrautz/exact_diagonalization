function H_up = calcHamiltonUpPointGroup( basisReprUp, compDownStatesPerRepr, ...
            paramT, indNeighbors, normHubbardStates, symOpInvariantsUp, kValue, ...
            index2ReprUp, symOp2ReprUp, intStatesUp, latticeSize)
%CALCHAMILTONIPPOINTGROUP %CALCHAMILTONUPONLYTRANS calculates the up spin 
% hopping part of the hubbard hamiltonian for only translational symmetry 
% and 1 and 2 dimensions
%
% Input:  basisReprUp             up spin representatives
%         compDownStatesPerRepr   compatible down spin states per up spin
%                                 representatives
%         paramT                  hopping parameter t
%         indNeigbors             index matrix of lattice neigbors.
%                                 lattice dimension implicetly included
%         normHubbardStates       norm of compatible hubbard basis states
%         symOpInvariantsUp       symmetry operations leaving 
%                                 representatives invariant: cell array!
%                                 cols 1:6 : symmetries with eigenvalues
%                                 col 7:     fermionic phase
%                                 col 8:     expontial phase
%         kValue                  k vector value
%         index2ReprUp            linking list of whole basis and UP spin
%                                 representatives. signed integer list! the
%                                 sign gives the phase factor of the state
%         symOp2ReprUp            symmetry operations linking basis states
%                                 corresponding representatives
%         intStatesUp             integer values of all reached Up spin
%                                 states. needed in calulation of up spin
%                                 hopping part for right linking in
%                                 combination with index2ReprUp
%
% Output: H_up                    up spin hopping part
%
% % NOTE: 23.05.13: changes from 1 dimensional instance of this function:
%   up until calculation of translationMatrix overlap code nearly exactly
%   the same. then have to extract X and Y translation component of symmetry
%   operations leaving UP spin invariant from symOpInvariantsUp and
%   symOp2ReprUp. change in respective indications and symmetry applications, and
%   calculation of exponential phases. plus have to apply new version of
%   translationMatrix overlap calculation.
% NOTE: 28.05.13: adaptation to new saving of down spin representatives
%  not much to do here, since i still save hubbard norm in old way. just
%  have to take care when i pick out compatible down spin parts
%  also made changes to memory efficiency
% NOTE: 03.06.13: found mistake which gave me wrong results for only
% translational symmetry comapered to no symmetry case: was in overlap
% matrix calculation <D'|Tx^jx Ty^jy|D> the power of translation was
% applied in wrong direction! changed in overlap cal. function to minus
% sign!
% NOTE: 04.06.13: change to point group:application of symOpInvariants! and
% symOp2ReprUp
% NOTE: 12.06.13: finally found all mistakes in this function!:
% 1) symmetry eigenvalues from mirroring were falsely included in
% cubicSymmetry function and so also in the phase in symOpInvariance so it
% was sometimes applied to often or in wrong places-> fixed!
% 2) also the function to calculate the overlap of DOWN spin states: the
% symmetry operations were in wrong order and symOp2ReprUp also in wrong
% place + here the symmetry eigenvalue of symOp2ReprUp also has to be
% applied!(wasnt before) is a problem because point group and translations
% do not commute all
% NOTE: 22.06.13: parallized second loop over connnected states and changed
% back to non referencing
%------------------------SVN Info------------------------------------------
% $Rev:: 82                                     $: Revision of last commit
% $Author:: dobrautz                            $: Author of last commit
% $Date:: 2013-06-22 19:51:27 +0200 (Sam, 22. J#$: Date of last commit
% -------------------------------------------------------------------------

% number of UP spin representatives and sites(1 and 2 dim)
[nReprUp, nSites] = size(basisReprUp);
% for binary2dez conversion
bin2dez = nSites-1:-1:0;
bin2dez = (2.^bin2dez)';
% determine if more than trivial T0 translation is possible
[nSymmetryOperations,~] = cellfun(@size, symOpInvariantsUp);
% cumulative index for Up spin Representatives, to calc. combined hubbard
% indices : UpIndex*numel(downStates) + DownIndex
% NOTE: 28.05.13: old form unfortunately had to change since new reference
% saving method, have to count now in each step
%!!! reversed it for norm of hubbard states! think of better way
cumulIndex = [0;cumsum(cellfun(@length,normHubbardStates))];
% number of dimensions of lattice:
[nDims,~] = size(indNeighbors);
% Translational phases to UP spin representatives
transPhases2ReprUp = sign(index2ReprUp);
% pull out phases of repr. linking list:
index2ReprUp = abs(index2ReprUp);

% container for matrix indices
% cell structure because of 2 dimensions
[xIndOfReprUp,yIndOfCycleUp,hoppingPhase] = deal(cell(nDims,1));

% -------------------------NOTE: 19.05.13----------------------------------
% this below is my known hopping function adapted for this function
%--------------------------------------------------------------------------

for iDims = 1:nDims %loop over Dimensions, her just one
    
    B2 = basisReprUp(:,indNeighbors(iDims,:)); % binary form of neighbors in DIM i
    B2 = xor(basisReprUp,B2);              % determines sites where electrons can hop
    % gets indices of sites where electrons can hop, but:
    [f,d] = find(B2);   %f zeilen, g spalten ,nach spalten sortiert oO
    
    % for UP Spin part one has to check if there is a repr. at all! or a repr.
    % with no hopping possibilities in a specific dimension
    if ~isempty(f)      %if there is hopping possible...
        
        % and a stupid matlab inconsitency has to be checked:
        % wenn nur eine zeile in B_repr macht mir find() einen zeilenvektor!!! wtf
        if nReprUp == 1
            f = f'; d = d'; % make column vector of it
        end
        
        [xIndOfReprUp{iDims},f] = sort(f);     % x_ind waere zeilen index von H_heisenberg
        
        % sort d same way as f, is first site index of 2 sites where hopping occurs
        d = d(f);
        f = indNeighbors(iDims,d)';          %i ndex of neighbour of flipped site
        
        f1 = [d;f];                       % combined index for flipped sites
        d1 = [1:sum(B2(:)),1:sum(B2(:))]';% two times index from 1:N for indexing the flipped sites
        % since i create in next row a new matrix
        % of binary form of, where basis is
        % reproduced so many times as hops can
        % occur
        
        F = basisReprUp(sum(B2,2)==0,:)*bin2dez;  %fehlende zahlen-indizes: as some
        % hops to integer values are not allowed
        % but one still needs this values for
        % correct ordering of new states to
        % determine y_ind value
        
        B2 = basisReprUp(xIndOfReprUp{iDims},:);  %basis so oft reproduziert wie gswitcht wird
        % NOTE: 23.05.13: need dummy variable here sometimes since it makes
        % indication matrix too small if there are zeros in wrong places
        d1 = sparse([d1;max(d1)],[f1;nSites],[ones(numel(d1),1);0]);         %indication matrix for switches
        
        % now calculate the fermionic phases of hamilton hops:
        % calc. hopping phases, maybe have to sort this if x and y are sorted
        [hoppingPhase{iDims}] = HubbardPhase(B2,d1,d,f);
        % TODO hubbardphase is bottleneck!! make it better!
        
        d = xor(d1,B2);                  % neue Zustaende!
        
        d = d*bin2dez;   % corresponding integer values
        f = [d;F];   % include missing integer values for comparison
        
        % also have to include integer values of cycle members, if they are
        % no reached by hoppings, needed for corrected linking with
        % index2ReprUp
        F = setdiff(intStatesUp,f);
        f = [f;F];%#ok
        
        [~,~,B2] = unique(f);% check them! thats enough!since standard call of unique
        % sorts the unique-values, indices B2 give
        % me right positions of new created states
        
        yIndOfCycleUp{iDims} = B2(1:numel(xIndOfReprUp{iDims}));%only include the one you need, possible since missing integer
        % values were included at the end of list
        
        
        clear f B2 d d1 f1
        
    end
end

% Indizes der nichtdiagonalelemente
xIndOfReprUp =cell2mat(xIndOfReprUp(:));
yIndOfCycleUp = cell2mat(yIndOfCycleUp(:));
hoppingPhase = cell2mat(hoppingPhase(:));
% transform H|U> to |U> representative:
% transform y index to representative
yIndOfReprUp = index2ReprUp(yIndOfCycleUp);

symOp2ReprUp = symOp2ReprUp(yIndOfCycleUp,:);
% combination of hopping phases and phase from translation to Repr.
combPhases = hoppingPhase .* transPhases2ReprUp(yIndOfCycleUp);

% transform indices to hubbard form, this + down spin part equals final
% index: UpIndex*numel(DownStates) + DownIndex
xIndHubbUp = cumulIndex(xIndOfReprUp);
yIndHubbUp = cumulIndex(yIndOfReprUp);

% bnumber of connected up spin states
nConnectedUpStates = numel(xIndOfReprUp);

% container for final x and y indices and phases
[xFinal,yFinal,phasesFinal] = deal(cell(nConnectedUpStates,1));

%-------------------------------17.05.13----------------------------------
% note: through nTranslations i know the number of indices to be spawned by
% a up spin state. maybe vectorization possible!
% NOTE: 19.05.13 not exactly because not only nTranslations but also the
% overlap of the translation matrix determines the number of spawned states
% -------------------------------------------------------------------------

% now loop over connected UP spin states( multiple sub indices possible)
% here changes in 2 dim instance of this
parfor iStates = 1:nConnectedUpStates;
    
    % hopped state index ( index of to be translated UP spin representative
    stateIndex = yIndOfReprUp(iStates);
    
    % pick out to compared DOWN spin parts, change from referencing:
    downSpinState1 = compDownStatesPerRepr{xIndOfReprUp(iStates)};%#ok
    downSpinState2 = compDownStatesPerRepr{stateIndex};
    
    % extract specific symmetry operations from symOpInvariant:
    % pick out symmetry ops which leave repr. invariant: 7th col:phase
    specSymmetryOperation = symOpInvariantsUp{stateIndex}(:,1:6);%#ok
    % phase from invariance sym op:
    phasesFromInvariance = symOpInvariantsUp{stateIndex}(:,7);
    
    % also have to take care of symmetry operation which bring hopped
    % states to correspondind representatives: take care of commutation
    % relations:!! TODO check which order they have to be applied!
    % go back to distinct function for PG and only translation case, is
    % more efficient
    [xIndHubbDown, yIndHubbDown, phaseFac] =  symOpMatrixOverlapPG( ...
        downSpinState1, downSpinState2,  specSymmetryOperation, ...
        symOp2ReprUp(iStates,:), latticeSize);
    
    % idea sketch
    xFinal{iStates} = xIndHubbUp(iStates) + xIndHubbDown;
    yFinal{iStates} = yIndHubbUp(iStates) + yIndHubbDown;
    
    % hopping and translation2Repr phase for spec. state
    % NOTE: 12.06.13:  combine it here with symmetry eigenvalue of
    % symOp2ReprUp: 
    symEigenValues2ReprUp  = symOp2ReprUp(iStates,3:end);
    symEigenValues2ReprUp = prod(symEigenValues2ReprUp(symEigenValues2ReprUp~=0));
    
    specCombPhase = combPhases(iStates)*symEigenValues2ReprUp;
    
    for iTrans = 1:nSymmetryOperations(stateIndex)%#ok
        
        % collect different phases
        phaseFromTransUp = phasesFromInvariance(iTrans);
        
        % symmetry eigenvalue of invariance
        symEigenValues = specSymmetryOperation(iTrans,3:end);
        symEigenValues = prod(symEigenValues(symEigenValues~=0));
    
        
        expPhases = exp(1i*(kValue(1)*(symOp2ReprUp(iStates,1) + specSymmetryOperation(iTrans,1)) + ...
                kValue(2)*(symOp2ReprUp(iStates,2) + specSymmetryOperation(iTrans,2)) ));%#ok
    
        % surpress real and imaginary parts of this below treshhold
        expPhases = real(expPhases)*(abs(real(expPhases))>10^-5) + ...
            1i*imag(expPhases)*(abs(imag(expPhases))>10^-5);
        
        % combine them
        phaseFac{iTrans} = phaseFac{iTrans}*specCombPhase*phaseFromTransUp ...
            *expPhases*symEigenValues;
        
    end
    % collect them
    phasesFinal{iStates} = cell2mat(phaseFac);
    
end

clearvars -EXCEPT *Final normHubbardStates paramT nSites

% convert cell arrays to vectors
maskIsEmpty = cellfun(@isempty,xFinal);
xFinal(maskIsEmpty) = [];
xFinal = cell2mat(xFinal); 
if ~isempty(xFinal)
    
    yFinal = cell2mat(yFinal);
    
    phasesFinal = cell2mat(phasesFinal);
    normHubbardStates = cell2mat(normHubbardStates);
    % number of states
    nHubbardStates = numel(normHubbardStates);
    % combine respective norms of states
    normHubbardStates = sqrt(normHubbardStates(xFinal) .* normHubbardStates(yFinal));
    % calculate hamiltonian elements
    H_up = -paramT .* phasesFinal ./ normHubbardStates ./ nSites;
    
    clear normHubbardStates phasesFinal
    
    % surpress values below treshhold
    H_up(abs(H_up)<10^-5) = 0;
    % include pseudo entry to ensure right size of hamiltonian
    [xFinal(end+1),yFinal(end+1)] = deal(nHubbardStates);
    % combine it in sparse form
    H_up = sparse(xFinal,yFinal,[H_up;0]);
    
    clear *Final

    % assure hermiticity
    H_up = (H_up + H_up')/2;

else
    
    H_up = sparse(1,1,0);
    
end


