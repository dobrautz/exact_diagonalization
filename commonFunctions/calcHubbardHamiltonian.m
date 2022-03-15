function H = calcHubbardHamiltonian(basisReprUp, compDownStatesPerRepr, ...
    compInd2ReprDown,paramU, paramT, indNeighbors, normHubbardStates,   ...
    symOpInvariantsUp, kValue, index2ReprUp, symOp2ReprUp, intStatesUp, ...
    basisStatesDown,latticeSize)
%POINTGROUPCALCHAMILTON calculation and combination of diagonal and UP/DOWN 
% spin hopping part of hubbard hamiltonian for all translations and one and 
% two dimensions
% 
% Input:    basisReprUp         Representatives of UP spin channel
%           compDownStatesPerRepr Representatives of DOWN spin channel in
%                               cell format. numel of this is the size of
%                               the hubbard basis. each cell element
%                               (matrix) corresponds to an UP spin repr.
%           paramU/T            parameters U and t of hubbard hamilton
%           indNeigbors         index matrix of lattice neighbors. this
%                               indicates dimension and lattice
%                               identification, for later generalization
%           normHubbardStates   cell array of norms of hubbard states
%           symOpInvariantsUp   symmetry operations leaving the UP spin
%                               part of the representatives invariant. also
%                               serves as identification which symmetries
%                               are applied!
%           kValue              k vector value
%           index2ReprUp        linking list of whole basis and UP spin
%                               representatives. signed integer list! the
%                               sign gives the phase factor of the state
%           symOp2ReprUp        symmetry operations linking basis states
%                               corresponding representatives
%           intStatesUp         integer values of all reached Up spin
%                               states. needed in calulation of up spin
%                               hopping part for right linking in
%                               combination with index2ReprUp
%           basisStatesDown     whole down spin basis
% 
% Output:   H                   Hubbard Hamilton Matrix
% 
% OLD NOTES:
% diagonal part of hamiltonian: <D|<U| Pk Hd |U>|D> stays diagonal, but in
% calculation one has to take care of reapearing states in Pk|Psi> and
% count how often thois happens for up anddown spin part and sum over resp.
% off-diagonal( hopping part) hamiltonian starts same way as i thought
% before but also have to apply projector Pk after Hh |psi> and than one
% has to look up ( calculate which states <U'|<D'| Pk Hh |U>|D> are
% connected: so one has to include the projector Pk in algorithms. or maybe
% there is a clever way of hash links to do the trick also! think about it!
% NOTE 13.05.2013:
% ad diagonal part of hamiltonian Hd:
% <rk|Hd|rk> = <r|Pk Hd Pk |r> / <r|Pk|r> = Hd(r) <r|Pk|r>/<r|Pk|r> = Hd(r)
% = summe of doubly occupancies!!! wait
% actually Hd(r) * <r|Pk|r>/ |<r|Pk|r>| = sgn(<Pk>)*Hd(r) is sign possible?
% TODO: check if sign possible!
% if sign is possible still easy to include. just use sign of normHubbard
% in right way in combination below!
% for now disregard possibility of negative or complex norm and go on:
% NOTE: 27.05.13: change in algorithm: to improve performance in down spin
% hopping matrix calcualtion: apply H|D> once to whole down spin basis and
% then link it to specific down spin part of representatives through
% compInd2ReprDown. can use this probably for other calculaions too:
% eg. make the transition overlap matrices for whole down spin basis and
% link to it through this list TODO
% NOTE: 28.05.13: have to adapt to new reference saving way of downState
% representatives and corresponding list
%  also made changes to memory efficiency
% NOTE: 04.06.13: inclusion of point group calculation to check if values
% of energies are right although number of basis states isnt correct yet!
% NOTE: 22.06.13: parallized a few code blocks in subfunctions from here
% also improved memory usage and hermiticity check also in subfunctions!
%------------------------SVN Info------------------------------------------
% $Rev:: 82                                     $: Revision of last commit 
% $Author:: dobrautz                            $: Author of last commit   
% $Date:: 2013-06-22 19:51:27 +0200 (Sam, 22. J#$: Date of last commit     
% -------------------------------------------------------------------------

% number of UP spin Repr 
[nReprUp,~] = size(basisReprUp);
nCompDownStates = sum(cellfun(@length,normHubbardStates));

if nReprUp == 0  || nCompDownStates == 0 % only do smth if there are representatives
    
    H = [];

else
    
    %% ------------------- Diagonal Term ----------------------------------
    % just count double occupancies!
    H = calcHubbardDIAG( basisReprUp, normHubbardStates, ...
        compDownStatesPerRepr, paramU);
    
    %% ----------------------- Off-Diagonal Term --------------------------
    %% ------------------------ DOWN Spin part ----------------------------
    % ------------------------- note 13.05.13------------------------------
    % <rk'|Hd(own)|rk> = N<r'|Hd|r> = N/L*sum(j=0..L-1) exp(2pi ijk/L) * 
    % <U'|Tj|U><D'|Tj Hd|D> -> diagonal in UP spin part!! <U|Tj|U> !!
    % Down spin part is diagonal in up spin part. since <U'|Tj|U> was
    % condition for chosing of representatives( two states which are only
    % different by a translation belong to same circle!) so further
    % calculation of down spin part is similar to my previous versions:
    % but not same amount of down spin states per repr. now!! but still
    % ordered correctly and everything. but for calculation i now also need
    % the connection <D'| Tj Hd |D> = hd(D)*<D'|Tj|D''> -> matrix Tj for
    % states needed! save it in previous function
    % plan:
    % loop over repr. ; application Hd|D>; loop over j: Tj|U> = pj|U> 
    % connection <D|Tj|D'> 
    %----------------------------------------------------------------------
    
    if ~iscell(symOpInvariantsUp)
    
        H = H + calcHamiltonDownOnlyTrans(compDownStatesPerRepr, compInd2ReprDown,...
            paramT, indNeighbors, normHubbardStates, symOpInvariantsUp, kValue,...
            basisStatesDown,latticeSize);
    
    else
        
        H = H + calcHamiltonDownPointGroup(compDownStatesPerRepr, compInd2ReprDown,...
            paramT, indNeighbors, normHubbardStates, symOpInvariantsUp, kValue,...
            basisStatesDown, latticeSize);
    end
    
    
    %% -------------------------- UP Spin Part ----------------------------
    % ----------------------------note 17.05.13----------------------------
    % <rk'|Hu(p)|rk> = N/L sum(j=0..L-1) e^(2pi ijk/L) <D'|Tj|D> *
    % <U'|Tj Hu |U> : now different calculation. because first of all Tj 
    % not diagonal in Down spin part <D'|Tj|D> ~= <D|Tj|D> and list is
    % sorted by UP spin part as MSB part. IDEA of calculation:
    % 1) calc. connection <U'|Tj Hu|U> = hu(U)<U'|Tj|U''> is kind of 
    % diagonal again. save tranlations connecting U'' and U' with indexRepr
    % list calculated earlier. and i have the list diagTransMatrixUP 
    % leaving U unchanged. save this information for UP spin state 
    % connections and phase factors indices etc. also
    % 2) loop over UP spin repr. and again over translations Tj within them
    % from point 1) and find the connection with the "small" Down spin
    % indices <D'|Tj|D>, here one has to take care of possible non
    % compatible states. some down spin states may be not common to both up
    % spin lists -> include in fct <D'|Tj|D>. with right combination this 
    % should be all
    %----------------------------------------------------------------------
    
    if ~iscell(symOpInvariantsUp)
        
        H = H + calcHamiltonUpOnlyTrans( basisReprUp, compDownStatesPerRepr, ...
            paramT, indNeighbors, normHubbardStates, symOpInvariantsUp, kValue, ...
            index2ReprUp, symOp2ReprUp, intStatesUp , latticeSize);
    else
        
        H = H + calcHamiltonUpPointGroup( basisReprUp, compDownStatesPerRepr, ...
            paramT, indNeighbors, normHubbardStates, symOpInvariantsUp, kValue, ...
            index2ReprUp, symOp2ReprUp, intStatesUp , latticeSize);
    end


end

