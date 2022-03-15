function diagHam = hamiltonDiagNoSyms( basisStatesUp, basisStatesDown,...
    paramU, flagSparse )
%HAMILTONDIAGNOSYMS calculates diagonal matrix elements of hubbard
%hamiltonian in binary repr. basis without any use of symmetry
% 
% Input: basisStatesUp/Down     binary repr. of heisenberglike basis
%        flagSparse             control flag if output should be in sparse
% 
% Output:diagHam                depending on flagSparse vector or sparse
%                               matrix output of diagonal matrix elements
% 
% NOTE: 07.10.13: created. to use it also in new FCIQMC algorithm

% number of UP, DOWN and Combined states
[nBasisStatesUp,~] = size(basisStatesUp);     
[nBasisStatesDown,~] = size(basisStatesDown);        

%% Diagonal Terms
% count double occupancies
diagHam = zeros(nBasisStatesDown,nBasisStatesUp); 

for iRepr = 1:nBasisStatesUp  
     % check which sites are doubly occupied with bsxfun not repmat
     isDoublyOccupied = bsxfun(@and,basisStatesUp(iRepr,:),basisStatesDown);
%      % replicate B_repr to size of B
%      repbasisStatesUp = repmat(basisStatesUp(iRepr,:),nbasisStatesDown,1);
%      % check which sites are doubly occupiedeyboard
%      isDoublyOccupied = and(repbasisStatesUp,basisStatesDown);
%      % sum over them. is this right? TODO
     diagHam(:,iRepr) = sum(isDoublyOccupied,2);
    
end
%diagonal terms
diagHam = paramU * diagHam(:);

if flagSparse == 1
    
    nHubbardStates = nBasisStatesUp*nBasisStatesDown; 
    
    diagHam =  sparse(1:nHubbardStates,1:nHubbardStates,diagHam);

end

