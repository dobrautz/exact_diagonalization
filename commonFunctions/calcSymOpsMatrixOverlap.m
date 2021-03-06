function [xInd, yInd, phaseFactor] = calcSymOpsMatrixOverlap( ...
    basisOne, basisTwo, flagSymmetry ,varargin)
%CALCSYMOPSMATRIXOVERLAP calculates overlap <D'|Tj Q|D> of symmetry
% combination defined by variable input elements varargin
% matrix of two basis sets D' and D for different values of j the two basis
% sets can have different sizes elements and even no overlap at all.
%
% Input:    basisOne/Two        the two to compared basis sets
%           flagSymmetry        indicating which case it is
%           (variable)          depending which symmetry different
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
% NOTE: 12.06.13: want to change this to a generla overlap function not 
% such a specified start as before! 
%------------------------SVN Info------------------------------------------
% $Rev:: 63                                     $: Revision of last commit
% $Author:: dobrautz                            $: Autor of last commit
% $Date:: 2013-06-13 09:05:27 +0200 (Don, 13. J#$: Datum of last commit
% -------------------------------------------------------------------------

% check for variable input variables
if strcmp(flagSymmetry,'1D')
    
    translationPower = varargin{1};
    % number of translation oerations
    nTranslations = numel(varargin{1});
    
elseif strcmp(flagSymmetry,'2D')
    
    transPowerX = varargin{1};
    transPowerY = varargin{2};
    % number of translation oerations
    nTranslations = numel(varargin{1});
    
elseif strcmp(flagSymmetry,'PG')
    
    specSymmetryOperation = varargin{1};
    symOp2ReprUp = varargin{2};
    
%     
%     specSymmetryOperation(:,1:2) = bsxfun(@plus,specSymmetryOperation(:,1:2),symOp2ReprUp(1:2));
%     specSymmetryOperation(:,3:end) = specSymmetryOperation(:,3:end).*...
%         bsxfun(@xor,specSymmetryOperation(:,3:end),symOp2ReprUp(3:end));
    
    % apply right inverse: negate translations:
    specSymmetryOperation(:,1:2) = -specSymmetryOperation(:,1:2);
    symOp2ReprUp(1:2) = -symOp2ReprUp(1:2);
    
    % number of symmetry ops
    [nTranslations,~] = size(varargin{1});
end

% keyboard
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
    
    switch flagSymmetry
        
        case '1D'
            
            % translated states and corresponding phases:
            % sign of tranlationPower not certain yet. for same results as in
            % fehske weisse paper positive sign!
            % NOTE: 03.06.13: -----------------------------------------------------
            % but this just a convention thing for right overlaps have to use
            % inverse -> negative signs!!
            [ transStates, transPhases] = cubicSymmetries(basisTwo, -translationPower(iTrans), latticeSize, 'rx1D');
            
        case '2D'
            % translated states and corresponding phases:
            % change in 2 dim case: first X translation
            [ transStates, transPhasesX] = cubicSymmetries(basisTwo, -transPowerX(iTrans), latticeSize, 'rx');
            % then Y translation
            [ transStates, transPhasesY] = cubicSymmetries(transStates, -transPowerY(iTrans), latticeSize, 'ry');
            
            transPhases = transPhasesX.* transPhasesY;
            
        case 'PG'
            
            [transStates,transPhasesX] = inverseSymmetry(basisTwo,specSymmetryOperation(iTrans,:), latticeSize);
            [transStates,transPhasesY] = inverseSymmetry(transStates,symOp2ReprUp, latticeSize);

            transPhases = transPhasesX .* transPhasesY;
    end
    
    
    % integer values
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


