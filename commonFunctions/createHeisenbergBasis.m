function [ basisStates, integerBasis, indOnes] = createHeisenbergBasis( nStates , nSites , S_ges )
% createHeisenbergBasis: Funktion zur Erstellung der Basis des Spin 1/2 -
% Heisenbergmodells
% basisStates = createHeisenbergBasis(nStates , nSites , S_ges) calculates
% the binary representation of a heisenberg hamiltonian
% [basisStates,integerBasis ] = createHeisenbergBasis(nStates,nSites,S_ges)
% also outputs the corresponding integer values of the binary form
% 
% Input:     nStates         number of states
%            nSites          number of sites
%            S_ges           Ordnung Gesamtspin
%            spar            Sparsamer Modus ein/aus
%
% Output:    basisStates     Basis in binaerer Darstellung
%            integerBasis    Basis in dezimaler Darstellung
%            indOnes         index list of occupied states
% 
% NOTE: 16.07.13: changed output to also be able to get indices of occupied
% states.
%------------------------SVN Info------------------------------------------
% $Rev:: 62                             $: Revision der letzten Übertragung
% $Author:: dobrautz                    $: Autor der letzten Übertragung
% $Date:: 2013-06-12 15:46:11 +0200 (Mi#$: Datum der letzten Übertragung
% -------------------------------------------------------------------------
% 
% TODO : check if matlab JIT Accelerator is running here

%% Basiserstellung
basisStates = false(nStates,nSites);           %container fuer Basis
%start mit kleinstmoeglicher zahl
basisStates(1,end:-1:(end-S_ges+1)) = 1; 

iSites = 1:nSites;                     %zur indexbestimmung der 1er und 0er

indOnes = zeros(nStates,sum(basisStates(1,:),2));

for iStates = 1:nStates-1
    % schleife zur berechnung der basis (geht laut VDL noch schneller)
    % wahrscheinlich in dem man pro schritt mehrere zahlen erzeugt
    maskOnes = basisStates(iStates,:) == 1;
    indexOnes = iSites(maskOnes);     %index der 1er
    indOnes(iStates,:) = indexOnes;
    indexZeros = iSites(~maskOnes);   %index der 0er
    
    % right most zero index, with atleast one right of it
    rightMostZero = max(indexZeros(indexZeros<max(indexOnes)));
    % duplicate i-1 basis
    basisStates(iStates+1,:) = basisStates(iStates,:);
    % set right most zero to one
    basisStates(iStates+1,rightMostZero) = 1;
    % and all right of it to zero
    basisStates(iStates+1,rightMostZero+1:end) = 0;
    % how many ones have been deleted
    nDeletedOnes = sum(indexOnes>rightMostZero+1);
    %set ones right of right most zero, with ones right of it
    basisStates(iStates+1,end-nDeletedOnes+1:end) = 1;
    
end

indOnes(end,:) = iSites(basisStates(end,:)==1);
if nargout > 1
    % zur umrechnung von binaer auf dezimal, overwrites matlab built-in bi2de!!
    % but this is faster anyways
    bi2de = nSites-1:-1:0;
    bi2de = 2.^bi2de;
    bi2de = bi2de';
    
    integerBasis = basisStates*bi2de;   %umrechnung bin2dec
end


