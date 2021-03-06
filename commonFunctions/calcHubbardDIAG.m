function H_diag = calcHubbardDIAG( basisReprUp, normHubbardStates, ...
    compDownStatesPerRepr, paramU)
%CALCHUBBARDDIAG calculation of diagonal part of Hubbard hamiltonian for
% all symmetries and one and two dimensions
% essentially only count double occupancies!
%
% Input:    basisReprUp             Representatives of UP spin channel
%           normHubbardStates       cell array of norms of hubbard states
%           compDownStatesPerRepr   Representatives of DOWN spin channel in
%           paramU                  parameters U of hubbard hamilton
% 
% Output:   H_diag                  diagonal part of hamiltonian
% 
% OLD NOTES:
% diagonal part of hamiltonian: <D|<U| Pk Hd |U>|D> stays diagonal 
% NOTE 13.05.2013:
% ad diagonal part of hamiltonian Hd:
% <rk|Hd|rk> = <r|Pk Hd Pk |r> / <r|Pk|r> = Hd(r) <r|Pk|r>/<r|Pk|r> = Hd(r)
% = summe of doubly occupancies!!! wait
% actually Hd(r) * <r|Pk|r>/ |<r|Pk|r>| = sgn(<Pk>)*Hd(r) is sign possible?
% TODO: check if sign possible!
% if sign is possible still easy to include. just use sign of normHubbard
% in right way in combination below!
% for now disregard possibility of negative or complex norm and go on:
% NOTE: 22.06.13: removed old check if states are only referencing and 
% parallized code
%------------------------SVN Info------------------------------------------
% $Rev:: 82                                     $: Revision of last commit 
% $Author:: dobrautz                            $: Author of last commit   
% $Date:: 2013-06-22 19:51:27 +0200 (Sam, 22. J#$: Date of last commit     
% -------------------------------------------------------------------------

% number of UP spin Repr and size of sepcific hubbard basis
[nReprUp,~] = size(basisReprUp);
nHubbardStates = sum(cellfun(@length,normHubbardStates));

% container for double occupancies
doubleOccupancies = cell(nReprUp,1);
%loop over repr.
parfor iRepr = 1:nReprUp
    
    % pick out down states per up spin repr.
    specBasisStatesDown = compDownStatesPerRepr{iRepr};
    
    % check which sites are doubly occupied with bsxfun not repmat
    isDoublyOccupied = bsxfun(@and, basisReprUp(iRepr,:), specBasisStatesDown);
    doubleOccupancies{iRepr} = sum(isDoublyOccupied,2);
    
end
% convert cell to vector
doubleOccupancies = cell2mat(doubleOccupancies);
matIndex = 1:nHubbardStates;

% diagonal terms <rk|Hd|rk> = <r|Pk Hd Pk |r> / <r|Pk|r> =
% = Hd(r) <r|Pk|r>/<r|Pk|r> = Hd(r);  see note above
H_diag = paramU * doubleOccupancies(:);
H_diag = sparse(matIndex,matIndex,H_diag);



