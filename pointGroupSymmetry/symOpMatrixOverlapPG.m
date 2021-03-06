function [xInd, yInd, phaseFactor] = symOpMatrixOverlapPG( ...
    basisOne, basisTwo, symOpOne, symOpTwo, latticeSize)
% symOpMatrixOverlapPG calculates overlap <D'|Tj Q|D> of symmetry
% operations including point groups and translations: the application is
% inversely to the left state !!! <D|(TjQ)^-1 !!
% matrix of two basis sets D' and D for different values of j the two basis
% sets can have different sizes elements and even no overlap at all.
%
% Input:    basisOne/Two        the two to compared basis sets
%           symOpOne/Two        symmetry operations between them
% Output:   xInd/yInd           indices connecting the two basis sets
%                               is a vector of size:(# connected states,1)
%           phaseFactors        phase factors associated with overlaps
%                               a cell array (# translationsPowers(
%                               # connected states,1),1)
%
% NOTE: 03.06.13: found mistake which gave me wrong results for only
% translational symmetry comapered to no symmetry case: was in overlap
% matrix calculation <D'|Tx^jx Ty^jy|D> the power of translation was
% applied in wrong direction! changed in overlap cal. function to minus
% sign!
% this change gives me opposite signs for imaginary parts as in Fehske
% Weisse paper! but right eigenvalues! probably just a convention thing
% NOTE: 04.06.13: implementation of point group overlap too, here one has
% to take care of commutation relations of operations!!
% NOTE: 12.06.13: distinct function for point group case! remember
% application is inversly to left state!! and also have to include symmetry
% eigenvalue of symOp2ReprUp(= symOpTwo for PG case!)
%------------------------SVN Info------------------------------------------
% $Rev:: 62                                     $: Revision of last commit
% $Author:: dobrautz                            $: Autor of last commit
% $Date:: 2013-06-12 15:46:11 +0200 (Mit, 12. J#$: Datum of last commit
% -------------------------------------------------------------------------


% apply right inverse: negate translations:
symOpOne(:,1:2) = -symOpOne(:,1:2);
symOpTwo(1:2) = -symOpTwo(1:2);
% number of symmetry ops
[nTranslations,~] = size(symOpOne);
% number of sites
[~,nSites] = size(basisOne);
% for binary2dez conversion
bin2dez = nSites-1:-1:0;
bin2dez = (2.^bin2dez)';
% integer Values of left states.
intStatesOne = basisOne*bin2dez;
% container variables for indces and phases
[xInd, yInd, phaseFactor] = deal(cell(nTranslations,1));

% loop over number of translations
for iTrans = 1:nTranslations
    
    % apply symmetry inversly and in inverse order too! to left state:
    [transStates,transPhasesX] = inverseSymmetry(basisTwo,symOpOne(iTrans,:), latticeSize);
    [transStates,transPhasesY] = inverseSymmetry(transStates,symOpTwo, latticeSize);
    
    transPhases = transPhasesX .* transPhasesY;

    intTransStates = transStates*bin2dez;
    
    % compare two integer lists: sketch
    % intersect() function does exactly this
    % and now non compatible representative indices dont even appear
    % since for intersect the values must be in both lists
    [~ , xInd{iTrans}, yInd{iTrans}] = intersect(intStatesOne,intTransStates);
    phaseFactor{iTrans} = transPhases(yInd{iTrans});
    
    
end

% convert to vector
xInd = cell2mat(xInd);
yInd = cell2mat(yInd);


