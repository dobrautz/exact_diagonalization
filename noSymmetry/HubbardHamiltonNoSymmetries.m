function H = HubbardHamiltonNoSymmetries( basisStatesUp, basisStatesDown, ...
    indNeighbors, paramU, paramT)
%HUBBARDHAMILTONNOSYMMETRIES calculates the full hubbard hamitlonian
% without application of symmetries for small system sizes
% NOTE: 16.06.13: tranperecy deteriorated a bit for better memory usage!
% renamed doubleOccupancies to H and change adding Up of H

% number of UP, DOWN and Combined states
[nBasisStatesUp,~] = size(basisStatesUp);     
      
%% Diagonal Terms
% count double occupancies
H = hamiltonDiagNoSyms( basisStatesUp, basisStatesDown, paramU, 1 );

%% Off-Diagonal Terms (seperate for UP and Down Spins)
%% DOWN (lsbs)
% (mit Algorithmus aus heisenberb bsp. : gibt mir indizes
% innerhalb der g representanten des up spin kanals
H = H + hamiltonDownNoSyms(basisStatesDown , nBasisStatesUp , paramT , indNeighbors);

%% UP (msbs)
%--------------------old notes---------------------------------------------
% Ordnung auch fix , da durch die basissortierung x = N*i + k alle UP repr.
% abstand N haben fuer fixes k ==> auch heisenberg anwendbar und dann index anpassung
%--------------------------------------------------------------------------
H = H + hamiltonUpNoSyms( basisStatesUp, basisStatesDown , paramT , indNeighbors );


