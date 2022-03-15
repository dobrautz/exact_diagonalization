function [downStatesPerRepr, index2ReprDown, normHubbardStates] = ...
    combine2HubbardBasisOnlyTrans(symOpInvariantsUp, basisStatesDown, latticeSize, kValue)
%COMBINE2HUBBARDBASIS combines Up and Down spin basis states to full
% hubbard basis for translational symmetry for one and two dimensions and
% calculates norm and compatibility of these states
%
% Input:        symOpInvariantsUp   symmetry operations leaving these
%                                   invariant
%               basisStatesDown     whole down spin basis set
%               nDims               number of dimension
%               kValue              value of curren k-vector
% 
% Output:       downStatesPerRepr   remaining down spin states per up spin
%                                   representative in cell format
%               index2ReprDown      list connecting each down states per up
%                                   spin repr. to whole down spin basis
%               normHubbardStates   norm of compatible(~=0) combined
%                                   Hubbard states
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
% NOTE: 12.06.13: this function not needed anymore! split functionalities
% again to improve performance and clarity
% NOTE: 21.06.13: removed referencing saving style of down spin parts or
% representatives to enable parallization
% NOTE: 22.06.13: included norm and compatibilty calculation of combined
% hubbard states also in this function like for point group symmetries for
% parallelization and cuntinuity reasons
% NOTE: 10.07.13: change latticeSize = [nSitesY;nSitesX] like matlab row
% column style
% 
%------------------------SVN Info------------------------------------------
% $Rev:: 82                                     $: Revision of last commit
% $Author:: dobrautz                            $: Author of last commit
% $Date:: 2013-06-22 19:51:27 +0200 (Sam, 22. J#$: Date of last commit
% -------------------------------------------------------------------------

% number of basis states of each spin component and number of sites
[nReprUp,~] = size(symOpInvariantsUp);
[nBasisStatesDown,nSites] = size(basisStatesDown);
[latticeDim,~] = size(latticeSize); 
nSitesY = latticeSize(1);

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

% loop over UP spin representatives
parfor iReprUp = 1:nReprUp
    
    % initialization of temporary parfor variables:
    transStates = 0;
    transPhases = 0;
    expPhase = 0;
    
    % for every UP spin representative determine invariance symmetries
    transIndexUp = find(symOpInvariantsUp(iReprUp,:));
    transPhasesUp = symOpInvariantsUp(iReprUp,transIndexUp);%#ok
    
    % most of the time there are no invariance symmetries except (0,0),
    % which means the whole down spin basis belongs to an UP spin repr.
    if numel(transIndexUp) == 1
        
        % for first time save it then reference to it
        downStatesPerRepr{iReprUp} = basisStatesDown;
        index2ReprDown{iReprUp} = indexDownStates;
        normHubbardStates{iReprUp} = ones(nBasisStatesDown,1)/nSites;
        
    else
        
        % dont need (0,0) element
        transIndexUp(1) = [];
        transPhasesUp(1) = [];
        
        % container reset for mask indicating that starting states are
        % smaller than symmetry applied states
        maskStatesSmaller = ones(nBasisStatesDown,1);
        % cummulatively add phases(symmetry,fermionic, exponential) up for
        % same staying states
        sumPhases = ones(nBasisStatesDown,1);
        
        % calculate combination of translational powers from:
        % columnPos-1 = translationPower; for one DIM
        translationPower = transIndexUp-1;
        % calculate combination of translational powers from:
        % columnPos = rx*Ly + ry + 1; for two DIM
        transPowerY = mod(transIndexUp-1,nSitesY);
        transPowerX = (transIndexUp - transPowerY - 1)/nSitesY;
        
        % otherwise apply translational symmetries, need loop because
        % cubicSymmetries cant handle vector symmetryvalue input;
        
        for iTrans = 1:numel(transIndexUp)
            
            if latticeDim == 1
                
                % apply translations: 1 dim translation:
                [transStates,transPhases] = cubicSymmetries(basisStatesDown,...
                    translationPower(iTrans), latticeSize, 'rx1D');
                
                expPhase = exp(1i*kValue*translationPower(iTrans)); 
                
            elseif latticeDim == 2
                
                % apply translations: first X translation:
                [transStates,transPhasesX] = cubicSymmetries(basisStatesDown,...
                    transPowerX(iTrans), latticeSize,'rx');
                % then Y translation
                [transStates,transPhasesY] = cubicSymmetries(transStates,...
                    transPowerY(iTrans), latticeSize,'ry');
                
                transPhases = transPhasesX .* transPhasesY;
                
                expPhase = exp(1i*(kValue(1)*transPowerX(iTrans) + ...
                    kValue(2)*transPowerY(iTrans)));
                
            end
            
            % integer values
            intTransStates = transStates*bin2dez;
            % update mask indicating if starting states are smaller than
            % symmetry applied states:
            maskStatesSmaller = and(maskStatesSmaller, (intDownStates <= intTransStates));
            
            % calc norm and compatibility:
            sameStates = intDownStates == intTransStates;
            
            sumPhases(sameStates) = sumPhases(sameStates) + expPhase * ...
                transPhasesUp(iTrans)*transPhases(sameStates);
            
        end
        
        % calc norm and compatibility:
        specNorm = abs(sumPhases)/nSites;
        maskStatesComp = specNorm > 10^-10;
        
        % combine logical compatibilty indicators( norm ~=0 and laFlorience
        % method of only keeping smaller states
        maskStatesComp = and(maskStatesComp, maskStatesSmaller);
        
        % only keep compatible and smaller states
        downStatesPerRepr{iReprUp} = basisStatesDown(maskStatesComp,:);
        index2ReprDown{iReprUp} = indexDownStates(maskStatesComp);
        normHubbardStates{iReprUp} = specNorm(maskStatesComp);
        
    end
end


