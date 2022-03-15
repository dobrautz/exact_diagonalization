function [exactDiagResults, H] = calcEigVectorsED(lattice, param, nSpins, exactDiagParam)
%CALCEIGVECTORSED function calculating ED eigenvectors for use in basis
%transformation function

% lattice size is minimal input: calc all from there since not costly

if ~isfield(lattice,'type')
    
    % assume cubic lattice as default
    lattice.type = 'cubic';
    
end

% parameters
% number of basis states for different spin channel
nStatesUp = nchoosek(lattice.nSites, nSpins.Up);
nStatesDown = nchoosek(lattice.nSites, nSpins.Down);

% Basis Creation
% without translational and point group symmetries
% NOTE: need integer values of whole UP spin basis for right linking in UP
% spin hopping part calculation of hamiltonian
[basisStatesUp, intStatesUp] = createHeisenbergBasis(nStatesUp, lattice.nSites, nSpins.Up);
[basisStatesDown, ~] = createHeisenbergBasis(nStatesDown, lattice.nSites, nSpins.Down);

if exactDiagParam.flagDoNoSymmetry && ~(lattice.nSites >=16 && (nSpins.Up == 8 && nSpins.Down == 8))
    
    tStartNoSym = tic;
    
    % calc hamiltonian
    H = HubbardHamiltonNoSymmetries(basisStatesUp, basisStatesDown, ...
        lattice.neighbors, param.U, param.t);
    
    % lanzcos
%     [exactDiagResults.noSymmetry.E0,exactDiagResults.noSymmetry.V0] = ...
%         lanczos_myversion(H,exactDiagParam);
%     
    [exactDiagResults.noSymmetry.V0,exactDiagResults.noSymmetry.E0] = ...
        eigs(H,1,'SA');
    
%     [a,b] = eig(full(H));
%     
%     exactDiagResults.noSymmetry.V0 = a(:,1);
%     exactDiagResults.noSymmetry.E0 = b(1);
%     
    exactDiagResults.noSymmetry.time = toc(tStartNoSym);

    exactDiagResults.noSymmetry.basisStatesUp = basisStatesUp;
    
    exactDiagResults.noSymmetry.basisStatesDown = basisStatesDown;
    
    exactDiagResults.noSymmetry.H = H;
%     clear H
end

if exactDiagParam.flagDoOnlyTrans && ~isempty(exactDiagParam.kNumber)
    
    tStartTrans = tic;
    % calculate representatives, symmetry operations leaving repr. invariant
    % and linking list and corresponding symmetry operations to repr.:
    [basisReprUp, symOpInvariantsUp, index2ReprUp, symOp2ReprUp ] = ...
        findReprOnlyTrans( basisStatesUp, lattice.size );
    
    % combine up and down spin states to hubbard basis laFlorience style
    % calculate the norm of hubbard basis states: here also incompatible
    % representatives per kValue are sorted out
    [compDownStatesPerRepr, compInd2ReprDown,normHubbardStates] = combine2HubbardBasisOnlyTrans(...
        symOpInvariantsUp, basisStatesDown, lattice.size, exactDiagParam.kValue);
    
    % calculate hamiltonian:
    H = calcHubbardHamiltonian(basisReprUp, compDownStatesPerRepr, ...
        compInd2ReprDown, param.U, param.t, lattice.neighbors, normHubbardStates,...
        symOpInvariantsUp, exactDiagParam.kValue, index2ReprUp, symOp2ReprUp, intStatesUp,...
        basisStatesDown,lattice.size);
    
    % lanczos
    [exactDiagResults.onlyTrans.E0, exactDiagResults.onlyTrans.V0] = ...
        lanczos_myversion(H,exactDiagParam);
    
    exactDiagResults.onlyTrans.norm = normHubbardStates;
    
    exactDiagResults.onlyTrans.time = toc(tStartTrans);
    
    exactDiagResults.onlyTrans.basisStatesUp = basisReprUp;
    
    exactDiagResults.onlyTrans.basisStatesDown = compDownStatesPerRepr;
    
    clear H comp* normHubbard* 
    
end

if exactDiagParam.flagApplyPointGroup && ~isempty(exactDiagParam.pointGroup) && lattice.dim == 2
    
    exactDiagParam.pointGroup = fliplr(exactDiagParam.pointGroup);

    tStartPoint = tic;
    % get right symmetries according do symmetry sector; CAUTION changes 
    % struct symmetrySector !!! 
    [whichSymmetry, exactDiagParam] = getSymmetryString(exactDiagParam, lattice);

    [basisReprUp, symOpInvariantsUp, index2ReprUp, symOp2ReprUp ] = ...
        PointGroupFindRepr( basisStatesUp, whichSymmetry, ...
        exactDiagParam.kValue, lattice.size, exactDiagParam.pointGroup);
    
    [compDownStatesPerRepr,compInd2ReprDown,normHubbardStates] = ...
        PointGroupLaFlorience( symOpInvariantsUp, basisStatesDown, lattice.size);
    
    H = calcHubbardHamiltonian(basisReprUp, compDownStatesPerRepr,...
        compInd2ReprDown,param.U, param.t, lattice.neighbors, normHubbardStates, ...
        symOpInvariantsUp,exactDiagParam.kValue, index2ReprUp, symOp2ReprUp, intStatesUp, ...
        basisStatesDown, lattice.size );
    
    
    [exactDiagResults.pointGroup.E0,exactDiagResults.pointGroup.V0] = ...
        lanczos_myversion(H,exactDiagParam);
    
    exactDiagResults.pointGroup.norm = normHubbardStates;
    
    exactDiagResults.pointGroup.time = toc(tStartPoint);
    
    exactDiagResults.pointGroup.basisStatesUp = basisReprUp;
    
    exactDiagResults.pointGroup.basisStatesDown = compDownStatesPerRepr;
    
end


    