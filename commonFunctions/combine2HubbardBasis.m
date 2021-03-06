function [ downStatesPerRepr, index2ReprDown, normHubbardStates ] = combine2HubbardBasis( ...
    symOpInvariantsUp, basisStatesDown,varargin)
%COMBINE2HUBBARDBASIS combines Up and Down spin basis states to full
% hubbard basis for translational and point group symmetry for one and two
% dimensions
%
% Input:        symOpInvariantsUp   symmetry operations leaving these
%                                   invariant
%               basisStatesDown     whole down spin basis set
%    (variable) nDims               number of dimension
%
% Output:       downStatesPerRepr   remaining down spin states per up spin
%                                   representative in cell format
%               index2ReprDown      list connecting each down states per up
%                                   spin repr. to whole down spin basis
%
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
% NOTE: 03.06.13: combination of translational and point group version of
% this function
% NOTE: 12.06.13: abandon this function and go back to distinct functions
% for each case! is faster and clearer
%------------------------SVN Info------------------------------------------
% $Rev:: 63                                     $: Revision of last commit
% $Author:: dobrautz                            $: Author of last commit
% $Date:: 2013-06-13 09:05:27 +0200 (Don, 13. J#$: Date of last commit
% -------------------------------------------------------------------------

if isempty(varargin)
    nDims = 2;
else
    nDims = varargin{1};
end

% keyboard
% number of basis states of each spin component and number of sites
[nReprUp,~] = size(symOpInvariantsUp);
[nBasisStatesDown,nSites] = size(basisStatesDown);
nSites1D = sqrt(nSites);
% umrechnung bin2dez
bin2dez = nSites-1:-1:0;
bin2dez = (2.^bin2dez)';
% integer Values of down spin basis States:
intDownStates = basisStatesDown*bin2dez;
indexDownStates = (1:nBasisStatesDown)';
% container for translated down spin states and linking list
downStatesPerRepr = cell(nReprUp,1);
index2ReprDown = cell(nReprUp,1);
normHubbardStates = cell(nReprUp,1);
% NOTE: 28.05.13: change way of saving index2ReprDown and downStatesPerRepr
% to save memory space because for 4x4 algorithm wont work:
flagAlreadySaved = 0;

% loop over UP spin representatives
for iReprUp = 1:nReprUp
    
    if iscell(symOpInvariantsUp)
        
        % NOTE: 24.05.13: change to PG here:
        % pick out specific symmetry operations: 7th column is phase factor!
        specSymOp = symOpInvariantsUp{iReprUp}(:,1:6);
        specSymPhase = symOpInvariantsUp{iReprUp}(:,7);
        specExpPhase = symOpInvariantsUp{iReprUp}(:,8);
        % number of sym. operations leaving it invariant (0,0) element always
        [nSpecSymOps,~] = size(specSymOp);
        
    else
        
        % for every UP spin representative determine invariance symmetries
        specSymOp = find(symOpInvariantsUp(iReprUp,:));
        specSymPhase = symOpInvariantsUp(iReprUp,specSymOp);
        nSpecSymOps = numel(specSymOp);
        
    end
    
    % most of the time there are no invariance symmetries except (0,0),
    % which means the whole down spin basis belongs to an UP spin repr.
    if nSpecSymOps == 1 && flagAlreadySaved == 0;
        
        % for first time save it then reference to it
        downStatesPerRepr{iReprUp} = basisStatesDown;
        index2ReprDown{iReprUp} = indexDownStates;
        normHubbardStates{iReprUp} = ones(nBasisStatesDown,1)/nSites;
        % change flag to position in cell:
        flagAlreadySaved = iReprUp;
        
    elseif nSpecSymOps == 1 && flagAlreadySaved ~= 0;
        
        % reference to cell position:
        downStatesPerRepr{iReprUp} = flagAlreadySaved;
        index2ReprDown{iReprUp} = flagAlreadySaved;
        normHubbardStates{iReprUp} = ones(nBasisStatesDown,1)/nSites;
    else
        
        % container reset for mask indicating that starting states are
        % smaller than symmetry applied states
        maskStatesSmaller = ones(nBasisStatesDown,1);

        if iscell(symOpInvariantsUp)
            
            % dont need (0,0) element
            specSymOp(1,:) = [];
            specSymPhase(1) = [];
            specExpPhase(1) = [];
            
            % cummulatively add phases(symmetry,fermionic, exponential) up for
            % same staying states
            sumPhases = ones(nBasisStatesDown,1);
            
            % apply symmetries: -1 because of cut of (0,0) element
            for iSymOp = 1:nSpecSymOps-1
                
                % apply spec symmetries:
                [symAppliedStates, symPhases] = symmetry(basisStatesDown, ...
                    specSymOp(iSymOp,:));
                
                % integer values
                intSymAppliedStates = symAppliedStates*bin2dez;
                
                % update mask indicating if starting states are smaller than
                % symmetry applied states:
                maskStatesSmaller = and(maskStatesSmaller, ...
                    (intDownStates <= intSymAppliedStates));
                
                % calc norm and compatibility:
                sameStates = intDownStates == intSymAppliedStates;
                
                symEigenValues = specSymOp(iSymOp,3:end);
                symEigenValues = prod(symEigenValues(symEigenValues~=0));
                
                
                sumPhases(sameStates) = sumPhases(sameStates) + specExpPhase(iSymOp) *...
                    symEigenValues* specSymPhase(iSymOp)* symPhases(sameStates);
                
            end
            
            specNorm = abs(sumPhases)/nSites;
            
            maskStatesComp = specNorm > 10^-10;
            
            maskStatesSmaller = and(maskStatesComp, maskStatesSmaller);
            
        else
            % dont need (0,0) element
            specSymOp(1) = [];
            
            % calculate combination of translational powers from:
            % columnPos-1 = translationPower; for one DIM
            translationPower = specSymOp-1;
            % calculate combination of translational powers from:
            % columnPos = rx*Ly + ry + 1; for two DIM
            transPowerY = mod(specSymOp-1,nSites1D);
            transPowerX = (specSymOp - transPowerY - 1)/nSites1D;
            
            % otherwise apply translational symmetries, need loop because
            % cubicSymmetries cant handle vector symmetryvalue input;
            for iSymOp = 1:nSpecSymOps-1
                
                if nDims == 1
                    
                    % apply translations: 1 dim translation:
                    transStates = cubicSymmetries(basisStatesDown,translationPower(iSymOp), latticeSize, 'rx1D');
                    
                elseif nDims == 2
                    
                    % apply translations: first X translation:
                    transStates = cubicSymmetries(basisStatesDown,transPowerX(iSymOp), latticeSize, 'rx');
                    % then Y translation
                    transStates = cubicSymmetries(transStates,transPowerY(iSymOp), latticeSize, 'ry');
                    
                end
                
                % integer values
                intTransStates = transStates*bin2dez;
                % update mask indicating if starting states are smaller than
                % symmetry applied states:
                maskStatesSmaller = and(maskStatesSmaller, ...
                    (intDownStates <= intTransStates));
                
            end
       
        end
        
        % only keep smaller states:
        downStatesPerRepr{iReprUp} = basisStatesDown(maskStatesSmaller,:);
        index2ReprDown{iReprUp} = indexDownStates(maskStatesSmaller);
        
        if ~iscell(symOpInvariantsUp)
            
            normHubbardStates{iReprUp} = specNorm(maskStatesSmaller);
            
        end
    end
end


