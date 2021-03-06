function [ symAppliedStates , symmetryPhases ] = inverseSymmetry( basisStates , symOperations , latticeSize)
%SYMMETRY applies symmetry operation s to binary basis and calcs
% fermionic phases
% Input:            basisStates ...    binary basis
%                   symOperations...   symmetry operation to apply
% 
% Output:           symAppliedStates   new symmetry applied basis states
%                   symmetryPhases     and corresponding fermionic phases
% 
% remember structure of sym_op: -rx -ry px py pd pe
% 
% group properties of 2D lattice-> can apply symmetries in a row:
% 
%------------------------SVN Info------------------------------------------
% $Rev:: 50                             $: Revision der letzten Übertragung
% $Author:: dobrautz                    $: Autor der letzten Übertragung
% $Date:: 2013-05-30 21:16:18 +0200 (Do#$: Datum der letzten Übertragung
% -------------------------------------------------------------------------
% 
% TODO: * bottleneck here just apply all symmetries since for eg px = 0, 
%           resp. fct does nothing actually 
symAppliedStates = basisStates;

%all symmetries: order of this operations not yet clear TODO
[symAppliedStates,phTy] = cubicSymmetries(symAppliedStates,symOperations(2),latticeSize,'ry');
[symAppliedStates,phTx] = cubicSymmetries(symAppliedStates,symOperations(1),latticeSize,'rx');
[symAppliedStates,phE] = cubicSymmetries(symAppliedStates,symOperations(6),latticeSize,'pe');
[symAppliedStates,phD] = cubicSymmetries(symAppliedStates,symOperations(5),latticeSize,'pd');
[symAppliedStates,phY] = cubicSymmetries(symAppliedStates,symOperations(4),latticeSize,'py');
[symAppliedStates,phX] = cubicSymmetries(symAppliedStates,symOperations(3),latticeSize,'px');


%--old version before inclusion of symmetry eigenvalue in phases phx from
% function cubicSymmetries.m
% symOperations = symOperations(3:end); 
% symOperations = symOperations(symOperations~=0);
% symmetryPhases = prod(symOperations) * phE.*phD.*phY.*phX.*phTx.*phTy;

symmetryPhases = phE.*phD.*phY.*phX.*phTx.*phTy;

end


