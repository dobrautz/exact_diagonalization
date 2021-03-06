function [downStatesPerRepr, index2ReprDown, normHubbardStates] =...
    PointGroupLaFlorience(symOpInvariantsUp, basisStatesDown, latticeSize)
%POINTGROUPLAFLORIENCE UP and Down Spin representative combination
% algorithm, according to LaFlorience method including point group
% symmetries, also calulates norm of hubbard states and checks
% compatibility
% 
% Input:        symOpInvariantsUp   matrix containing symmetry operations 
%                                   leaving each UP spin representatives 
%                                   invariant. contains symOps like: 
%                                   kx ky px py pd pe as eigenvalues and
%                                   fermionic and exp. phase in 7th, 8th col.
%               basisStatesDown     whole down spin basis set
%
% Output:       downStatesPerRepr   remaining down spin states per up spin
%                                   representative in cell format
%               index2ReprDown      list connecting each down states per up
%                                   spin repr. to whole down spin basis 
% 
% NOTE: 24.05.13: copied together from only translation version of this: 
%   difference is that symOpInvariantsUp is now a CELL and different
%   symmetries have to be applied to the basisStatesDown
% NOTE: 22.05.13: difference in LaFlorience method in basis combination:
% only calculate representatives for UP spin part and keep track of
% symmetry operations leaving UP spin representatives invariant. ( already
% have this list) than make tensorial product of |U_repr> x |D_all> than
% apply aforementioned symmetries to down spin part and only keep those
% with smaller integer value.
% NOTE: 28.05.13: change way of saving index2ReprDown and downStatesPerRepr
% to save memory space because for 4x4 algorithm wont work:
% only reference to a state if down spin part would be whole down spin
% basis: would be same information
% NOTE: 12.06.13: probably possible to calculate norm of hubbard basis in
% here too! since i need phases from same staying elements-> need down spin
% states which stay the same under up spin symmetry invariances! and sum
% over their phases to see if they are compatible! but doesnt work yet!
% wouldnt need hubbardBasisNorm function! DONE
%------------------------SVN Info------------------------------------------
% $Rev:: 82                                     $: Revision of last commit
% $Author:: dobrautz                            $: Author of last commit
% $Date:: 2013-06-22 19:51:27 +0200 (Sam, 22. J#$: Date of last commit
% -------------------------------------------------------------------------

% number of basis states of each spin component and number of sites
[nReprUp,~] = size(symOpInvariantsUp);
[nBasisStatesDown,nSites] = size(basisStatesDown);
% umrechnung bin2dez
bin2dez = nSites-1:-1:0;
bin2dez = (2.^bin2dez)';
% integer Values of down spin basis States:
intDownStates = basisStatesDown*bin2dez;
indexDownStates = (1:nBasisStatesDown)';
% container for translated down spin states
downStatesPerRepr = cell(nReprUp,1);
index2ReprDown = cell(nReprUp,1);
normHubbardStates = cell(nReprUp,1);

% loop over UP spin representatives
parfor iReprUp = 1:nReprUp
    
    % NOTE: 24.05.13: change to PG here:
    % pick out specific symmetry operations: 7th column is phase factor!
    specSymOp = symOpInvariantsUp{iReprUp}(:,1:6); 
    specSymPhase = symOpInvariantsUp{iReprUp}(:,7);
    specExpPhase = symOpInvariantsUp{iReprUp}(:,8);
    
    % number of sym. operations leaving it invariant (0,0) element always
    [nSpecSymOps,~] = size(specSymOp);
    
    % most of the time there are no invariance symmetries except (0,0),
    % which means the whole down spin basis belongs to an UP spin repr.
    if nSpecSymOps == 1
        
        downStatesPerRepr{iReprUp} = basisStatesDown;
        index2ReprDown{iReprUp} = indexDownStates;
        normHubbardStates{iReprUp} = ones(nBasisStatesDown,1)/nSites;
            
    else % apply corresponding symmetry operations
        
        % dont need (0,0) element
        specSymOp(1,:) = [];
        specSymPhase(1) = [];
        specExpPhase(1) = [];
        
        % specSymOp now contains: kx ky px py pd pe to be applied
        
        % container reset for mask indicating that starting states are
        % smaller than symmetry applied states
        maskStatesSmaller = ones(nBasisStatesDown,1);
        
        % cummulatively add phases(symmetry,fermionic, exponential) up for
        % same staying states
        sumPhases = ones(nBasisStatesDown,1);
        
        % apply symmetries: -1 because of cut of (0,0) element
        for iSymOp = 1:nSpecSymOps-1
            
            % apply spec symmetries:
            [symAppliedStates, symPhases] = symmetry(basisStatesDown, specSymOp(iSymOp,:), latticeSize);
            
            % integer values 
            intSymAppliedStates = symAppliedStates*bin2dez;
            
            % update mask indicating if starting states are smaller than
            % symmetry applied states:
            maskStatesSmaller = and(maskStatesSmaller, (intDownStates <= intSymAppliedStates));
            
            % calc norm and compatibility:
            sameStates = intDownStates == intSymAppliedStates;
            
            symEigenValues = specSymOp(iSymOp,3:end);
            symEigenValues = prod(symEigenValues(symEigenValues~=0));
            
            
            sumPhases(sameStates) = sumPhases(sameStates) + specExpPhase(iSymOp) *...
                symEigenValues* specSymPhase(iSymOp)* symPhases(sameStates);
            
        end
        
        specNorm = abs(sumPhases)/nSites;
        
        maskStatesComp = specNorm > 10^-10;
        
        maskStatesComp = and(maskStatesComp, maskStatesSmaller);
        
        % collect compatible states and hash list
        downStatesPerRepr{iReprUp} = basisStatesDown(maskStatesComp,:);
        index2ReprDown{iReprUp} = indexDownStates(maskStatesComp);
        normHubbardStates{iReprUp} = specNorm(maskStatesComp);
        
    end
end


