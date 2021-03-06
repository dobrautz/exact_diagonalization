function [ xInd, yInd, phaseFactor ] = calcTransMatrOverlapOnlyTrans( ...
    basisOne, basisTwo, latticeSize, varargin)
%CALCTRANSMATROVERLAPONLYTRANS calculates overlap <D'|Tj|D> of translation 
% matrix of two basis sets D' and D for different values of j the two basis
% sets can have different sizes elements and even no overlap at all.
% 
% Input:    basisOne/Two        the two to compared basis sets   
%(variable) translationPower    translational power j of T^j = how many 
%                               translations 
% 
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
% 
%------------------------SVN Info------------------------------------------
% $Rev:: 54                                     $: Revision of last commit
% $Author:: dobrautz                            $: Autor of last commit
% $Date:: 2013-06-03 11:26:59 +0200 (Mon, 03. J#$: Datum of last commit
% -------------------------------------------------------------------------

% check for variable input variables
if numel(varargin) == 1
    translationPower = varargin{1};
    nDims = 1;
    
elseif numel(varargin) == 2
    transPowerX = varargin{1};
    transPowerY = varargin{2};
    nDims = 2;
end

% number of sites
[~,nSites] = size(basisOne);
% number of translation oerations
nTranslations = numel(varargin{1});
% for binary2dez conversion
bin2dez = nSites-1:-1:0;
bin2dez = (2.^bin2dez)';
% integer Values of left states.
intStatesOne = basisOne*bin2dez;
% container variables for indces and phases
[xInd, yInd, phaseFactor] = deal(cell(nTranslations,1));

% loop over number of translations
for iTrans = 1:nTranslations
    
    if nDims == 1
        
    % translated states and corresponding phases:
    % sign of tranlationPower not certain yet. for same results as in
    % fehske weisse paper positive sign!
    % NOTE: 03.06.13: -----------------------------------------------------
    % but this just a convention thing for right overlaps have to use
    % inverse -> negative signs!!
    [ transStates, transPhases] = cubicSymmetries(basisTwo, -translationPower(iTrans), latticeSize, 'rx1D');
    
    elseif nDims == 2
    % translated states and corresponding phases:
    % change in 2 dim case: first X translation
    [ transStates, transPhasesX] = cubicSymmetries(basisTwo, -transPowerX(iTrans), latticeSize, 'rx');
    % then Y translation
    [ transStates, transPhasesY] = cubicSymmetries(transStates, -transPowerY(iTrans), latticeSize, 'ry');
    
    transPhases = transPhasesX.* transPhasesY;
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

