function H_down = calcHamiltonDownOnlyTrans(compDownStatesPerRepr, compInd2ReprDown,...
        paramT, indNeighbors, normHubbardStates, symOpInvariantsUp, kValue,...
        basisStatesDown,latticeSize)
%CALCHAMILTONDOWNONLYTRANS calculates the down spin hopping part of the
%hubbard hamiltonian for only translational symmetry and 1 and 2 dimensions
% 
% Input:    compDownStatesPerRepr   compatible down spin states per up spin
%                                   representatives
%           compInd2ReprDown        list linking down spin representative
%                                   to whole down spin basis
%           paramT                  hopping parameter t
%           indNeigbors             index matrix of lattice neigbors.
%                                   lattice dimension implicetly included
%           normHubbardStates       norm of compatible hubbard basis states
%           symOpInvariantsUp       matrix containing power and phase
%                                   factor leaving up spin representatives
%                                   invariant under translations. non zero
%                                   elements are phases columnPosition-1 is
%                                   translation power
%           kValue                  k vector value
%           basisiStatesDown        whole down spin basis
% 
% Output:   H_down                  Down spin hopping part
% 
% NOTE: 19.05.13 :improve code: maybe to remove loop over up spin 
% representatives calulate T^j H |D> once for whole heisenberg-like down 
% spin basis set and work with linking list to it for down spin parts
% of representatives; DONE!
% NOTE: 23.05.13: changes from 1 dimensional instance:
%       have to extract X and Y translation component of symmetry
%       operations leaving UP spin invariant from symOpInvariantsUp and
%       change respective indications and symmetry applications, and
%       calculation of exponential phases. rest identical
% NOTE: 23.05.13: apperently there are now diagonal elements in hopping
%   hamiltonians! check if this is right TODO!
% NOTE: 27.05.13: change in algorithm: to improve performance in down spin
% hopping matrix calcualtion: apply H|D> once to whole down spin basis and
% then link it to specific down spin part of representatives through
% compInd2ReprDown.
% NOTE: 28.05.13: would be better to change algorithm even more and also
% calculate translations on whole down spin basis and then sort out
% incompatible sates with linking lists, only problems now
% NOTE: 28.05.13: adaptation of new reference saving method: also huge
% improvement here: if specific states are whole basis -> there are no
% additional symmetries then (0,0) -> have just 'intern' spin down index
% distribution for each up spin -> as before(long ago)
% first part: H application to whole basis stays the same:
%  also made changes to memory efficiency
% NOTE: 18.06.13: changed ismember() call to undocumented function ismembc!
% NOTE: 22.06.13: parallized second loop over up representatives
%------------------------SVN Info------------------------------------------
% $Rev:: 82                                     $: Revision of last commit 
% $Author:: dobrautz                            $: Author of last commit   
% $Date:: 2013-06-22 19:51:27 +0200 (Sam, 22. J#$: Date of last commit     
% -------------------------------------------------------------------------

% number of UP spin representatives and sites( 1 and 2 Dim)
[nReprUp, nSites] = size(symOpInvariantsUp);
nSitesY = latticeSize(1);

% for binary2dez conversion
bin2dez = nSites-1:-1:0;
bin2dez = (2.^bin2dez)';
% determine if more than trivial T^0 translation is possible
nTranslations = sum(symOpInvariantsUp~=0,2);
% cumulative index for Up spin Representatives, to calc. combined hubbard
% indices : UpIndex*numel(downStates) + DownIndex
cumulIndex = [0;cumsum(cellfun(@length,normHubbardStates))];
% cumulIndex = 0;
% number of dimensions of lattice:
[nDims,~] = size(indNeighbors);
% number of basis states of whole down spin basis
[nBasisSatesDown,~] = size(basisStatesDown);
% final x and y indices and phases
[xFinal,yFinal,phasesFinal] = deal(cell(nReprUp,1));

%container for matrix indices (need cell structure because of 2 dimensions)
[xIndWholeBasis,yIndWholeBasis,phasesWholeBasis] = deal(cell(nDims,1)); 
%% application of hamiltonian to whole down spin basis

for iDims = 1:nDims %loop over Dimensions
    
B2 = basisStatesDown(:,indNeighbors(iDims,:));      %binary form of neighbors in DIM i
B2 = xor(basisStatesDown,B2);         % determines sites where electrons can hop

%gets indices of sites where electrons can hop, but:
[f,d] = find(B2);         %f zeilen, d spalten ,nach spalten sortiert oO
[xIndWholeBasis{iDims},f] = sort(f);   % x_ind waere zeilen index von H_heisenberg

d = d(f);                 %sort d same way as f, is first site index of 2 sites where hopping occurs
f = indNeighbors(iDims,d)';             %index of neighbour of to be flipped site
f1 = [d;f];               %combined index for flipped sites
d1 = [1:sum(B2(:)),1:sum(B2(:))]';%two times index from 1:N for indexing the flipped sites
                                  %since i create in next row a new matrix
                                  %of binary form of, where basis is
                                  %reproduced so many times as hops can
                                  %occur
                                  
F = basisStatesDown(sum(B2,2)==0,:)*bin2dez;           %fehlende zahlen-indizes: as some hops to integer values are not allowed
                                    %but one still needs this values for
                                    %correct ordering of new states to
                                    %determine y_ind value
                                    
B2 = basisStatesDown(xIndWholeBasis{iDims},:);     %basis so oft reproduziert wie gswitcht wird, see above

% NOTE: 23.05.13: need dummy variable here sometimes since it makes
% indication matrix too small if there are zeros in wrong places
d1 = sparse([d1;max(d1)],[f1;nSites],[ones(numel(d1),1);0]); %indication matrix for switches

% now calculate the fermionic phases of hamilton hops:
[phasesWholeBasis{iDims}] = HubbardPhase(B2,d1,d,f);%maybe have to sort this if x and y are sorted
%TO DO hubbardphase is bottleneck!! make it better!

d = xor(d1,B2);                  %new states, thats enough

d = d*bin2dez;                       %corresponding integer values
f = [d;F];                      %include missing integer values for comparison

[~,~,B2] = unique(f);           %check them! thats enough! since standard call of unique
                                % sorts the unique-values, indices B2 give
                                % me right positions of new created states 
yIndWholeBasis{iDims} = B2(1:numel(xIndWholeBasis{iDims}));       %only include the one you need, possible since missing integer
                                        % values were included at the end
                                        % of list
                         
clear f B2 d d1 f1 

end

%Indizes der nichtdiagonalelemente
xIndWholeBasis = cell2mat(xIndWholeBasis(:));  
yIndWholeBasis = cell2mat(yIndWholeBasis(:)); 
phasesWholeBasis = cell2mat(phasesWholeBasis(:)); 

% loop over UP spin representatives
parfor iReprUp = 1:nReprUp
    
    % intitialize temporary parfor variables
    transBasis = 0;
    transPhases = 0;
    expPhases = 0;
    
    % specific number of translations per UP spin repr.
    specNumOfTrans = nTranslations(iReprUp);
    
    % if there is only one translation = (0,0)-one : do nothing except
    % cumulatively add number of down spin states for indicating different
    % up spin states:
    if specNumOfTrans == 1 
        % put in indices and only convert them with corresponding 
        % cumulative index:
        xFinal{iReprUp} = cumulIndex(iReprUp) + xIndWholeBasis; 
        yFinal{iReprUp} = cumulIndex(iReprUp) + yIndWholeBasis; 
        phasesFinal{iReprUp} = phasesWholeBasis;
        
        % update cumulative Index
%         cumulIndex = cumulIndex + nBasisSatesDown;

    else
    
    % container for matrix indices
    % cell structure because of 2 dimensions and different numbers of
    % invariance translations
    [xIndOfReprDown,yIndOfReprDown,hoppingPhase] = deal(cell(specNumOfTrans,1));
    
    % pick out secific linking list between whole basis and specific one:
    specInd2ReprDown = compInd2ReprDown{iReprUp};
    
    % number of specific down states
    nSpecStatesDown = numel(specInd2ReprDown);

    % make index tranformation from whole basis indication to specific one:
    indexTransform = zeros(nBasisSatesDown,1);
    indexTransform(specInd2ReprDown) = 1:nSpecStatesDown;
    % transform it and exclude non compatible ones
    xIndSpec = indexTransform(xIndWholeBasis);
    yIndSpec = indexTransform(yIndWholeBasis);
    
    % mask specific possible starting states:
    maskStartingStates = xIndSpec~=0;
    
    % without translations: compatible end states with linking list:
    maskCompatible = and(maskStartingStates, yIndSpec ~= 0);
    
    % copy indices and phases before i cut it for (0,0) translation
    if specNumOfTrans > 1
        for iTrans = 2:specNumOfTrans
            
            xIndOfReprDown{iTrans} = xIndSpec(maskStartingStates);
            hoppingPhase{iTrans} = phasesWholeBasis(maskStartingStates);
        end
    end
    % cut it
    xIndSpec = xIndSpec(maskCompatible);
    yIndSpec = yIndSpec(maskCompatible);
    phasesSpec = phasesWholeBasis(maskCompatible);
   
     % j = 0 translation power is always: to this first
    xIndOfReprDown{1} = xIndSpec;
    yIndOfReprDown{1} = yIndSpec;
    hoppingPhase{1} = phasesSpec;
    % thats all if there is only (0,0) translation;
    
    % ------------------------ NOTE 13.05.13---------------------------
    % until here its the same as in the 2dim version of this function,
    % but now some differences are coming:
    % 2) here has the loop over translation powers j has to occur. TODO
    % find correct implication of Tj|d> and linking to starting state..
    %------------------------------------------------------------------
    
    if specNumOfTrans > 1
        
        % pick out Down spin basis states and calc. integers
        specStatesDown = compDownStatesPerRepr{iReprUp};
        intSpecStatesDown = specStatesDown*bin2dez;
        
        % recreate binary basis from within hamilton calculation:
        % NOTE: 28.05.13: have to take states from possible starting states
        % for specific Up spin Repr. after  hopping, can be states outside
        % of specific down states! 
        hoppedStates = basisStatesDown(yIndWholeBasis(maskStartingStates),:);%#ok
        
        % translational powers leavin up spin repr. invariant
        [~,translationPower,phasesFromTransInvariance] = find(symOpInvariantsUp(iReprUp,:));
        
        % calculate combination of translational powers from:
        % columnPos = rx*Ly + ry + 1; for two DIM
        transPowerY = mod(translationPower-1,nSitesY);
        transPowerX = (translationPower - transPowerY - 1)/nSitesY;
        % translational powers leavin up spin repr. invariant for 1 DIM
        translationPower = translationPower-1;
        
        for iTrans = 2:specNumOfTrans
            
            % pick out specific translations 2D
            specTransX = transPowerX(iTrans);
            specTransY = transPowerY(iTrans);
            
            % pick out specific translation 1D
            specTrans = translationPower(iTrans);
            
            if nDims == 1
                
                % apply translation
                [transBasis, transPhases] = cubicSymmetries(hoppedStates, specTrans , latticeSize, 'rx1D');
                % collect different phases
                expPhases = exp(1i*kValue*specTrans);
                
            else
                
                % apply translations: first X-translation:
                [transBasis, transPhasesX] = cubicSymmetries( hoppedStates, specTransX , latticeSize, 'rx');
                % then Y translation:
                [transBasis, transPhasesY] = cubicSymmetries( transBasis, specTransY , latticeSize, 'ry');
                % collect different phases
                expPhases = exp(1i*(kValue(1)*specTransX + kValue(2)*specTransY));
                transPhases = transPhasesX.*transPhasesY;
                
            end
            
            
            % surpress real and imaginary parts of this below treshhold
            expPhases = real(expPhases)*(abs(real(expPhases))>10^-5) + ...
                1i*imag(expPhases)*(abs(imag(expPhases))>10^-5);
            % phase from up spin repr. translational invariance
            phaseUpSpinTransInv = phasesFromTransInvariance(iTrans);
            
            % combine all obtained phase factors to one:
            hoppingPhase{iTrans} = hoppingPhase{iTrans}.*transPhases.* ...
                phaseUpSpinTransInv* expPhases;
            
            % sort out incompatible reached states: ismember() fct!
            intValues = transBasis*bin2dez;
            maskInBasis = ismembc(intValues,intSpecStatesDown);
%             maskInBasis = ismember(intValues,intSpecStatesDown);

            % only keep those which give rise to compatible representatives
            intValues = intValues(maskInBasis);
            xIndOfReprDown{iTrans}(~maskInBasis) = [];
            hoppingPhase{iTrans}(~maskInBasis) = [];
            
            % include missing integer values for comparison
            % now more careful here: through translation new states are
            % created. have to include states which are in Basis but not
            % in transBasis
            F = setdiff(intSpecStatesDown, intValues);
            f = [intValues;F];
            [~,~,B2] = unique(f);
            % TODO find better way to handle cases where no state is
            % reached
            yIndOfReprDown{iTrans} = B2(1:numel(xIndOfReprDown{iTrans}));
            
        end
    end
    
    % combine UP and DOWN spin index factors: 
    % UpInd*numel(downStates) + DownInd
    % NOTE: 23.05.13: have to tke care of cases where no states are
    % connected, compatible or otherwise: cut those out
    maskIsEmpty = cellfun(@isempty,xIndOfReprDown);
    xIndOfReprDown(maskIsEmpty) = [];
    yIndOfReprDown(maskIsEmpty) = [];
    % convert them with corresponding cumulative index:
    xFinal{iReprUp} = cumulIndex(iReprUp) + cell2mat(xIndOfReprDown(:)); 
    yFinal{iReprUp} = cumulIndex(iReprUp) + cell2mat(yIndOfReprDown(:)); 
    phasesFinal{iReprUp} = cell2mat(hoppingPhase(:));
    
    % update cumulative Index:
%     cumulIndex = cumulIndex + nSpecStatesDown;
    end
end

clearvars -EXCEPT *Final normHubbardStates paramT nSites

% convert cell arrays to vectors
xFinal = cell2mat(xFinal(:)); 
yFinal = cell2mat(yFinal(:));
phasesFinal = cell2mat(phasesFinal(:));
normHubbardStates = cell2mat(normHubbardStates);
% number of states
nHubbardStates = numel(normHubbardStates);
% combine respective norms of states
normHubbardStates = sqrt(normHubbardStates(xFinal) .* normHubbardStates(yFinal));

% calculate hamiltonian elements
H_down = -paramT / nSites .* phasesFinal ./ normHubbardStates;

clear normHubbardStates phasesFinal

% surpress values below treshhold
H_down(abs(H_down)<10^-5) = 0;

% combine it in sparse form
% include pseudo entry to ensure right size of hamiltonian
H_down = sparse([xFinal;nHubbardStates],[yFinal;nHubbardStates],[H_down;0]);

clear *Final

% assure hermiticity
H_down = (H_down + H_down')/2;

