%% Exact Diagonalization Algorithm for the Hubbard Model
% !!! Important normal order := first spin then lattice index !!!
% 
% for now: just works for 1 dim and only translation symmetry but can 
% include non half filling case.
% 
% first version:     10.12.2012
% second version:    09.04.2013 (major update -> norm of repr.)
% 1st branch:        10.05.2013 include non half filling and 1Dim
%                    18.05.13   final seems to work for 1 Dim
% 
%------------------------SVN Info------------------------------------------
% $Rev:: 87                                     $: Revision of last commit
% $Author:: dobrautz                            $: Author of last commit
% $Date:: 2013-06-25 10:04:48 +0200 (Die, 25. J#$: Date of last commit
% -------------------------------------------------------------------------
% 
%--------------------- info 12.04.2013:------------------------------------
% test run with Lx = 4, dauert  23510sec = 6,5h
% matrixgroesse fuer letzten block kx = ky = pi/2 war 4993560^2
% ueberpruefe ob ich eine so grosse matrix mit komplexen zahlen ueberhaupt
% aufstellen kann, da diese matrix 171020556 eintraege hat!
% TODO:
% noch einen test run machen und matrixgroesse und # eintrage aller machen
%--------------------------------------------------------------------------
% 
%----------------------------note 29.04.2013-------------------------------
% Get error:
% In an assignment  A(I) = B, the number of elements in B and I must be the same.
% Error in calcRepresentants (line 155)
%         indexRepr(iBasis+1) = alreadyRepr;  %position of corresponding repr.
% Error in applyCubicLatticeSymmetries (line 314)
% [ indexRepr , symOperationsRepr , basisRepr , normRepr] = ...
% Error in Hubbard_ED (line 63)
%             [indexRepr , symOperations2Repr , basisRepr  , normRepr] = ... 
% 
% but only if i run Hubbard_ED without breakpoints or keyboard keywords...
% if i run it step by step in debugger it finishes without an error!! oO
% 
%------------------update 30.04.13 on note above---------------------------
% found bug. seems like its a matlab bug. in the second and further runs of
% the for loops within the if statements, somehow matlab doesnt recognise 
% the inputname anymore in the function calcRepresentants.m ??? 
% TODO reproduce this error in a smaller function and send it to mathworks
%--------------------------------------------------------------------------
% 
% TODO: * reproduce matlab inputname() bug
%       * do test run with TAG release which works for Lx = 4 and save sizes
%           and number of eintraege of hamilton matrix, and beforehand update
%           the their included TableUpdate.m
%          + check hermiticity of the 4x4 lattice hamiltonians
%       * read papers:  ED Techniques; Weisse, Fehske
%                       Space Group Symmetries; Laflorencie Poilblanc 
%           concerning combination of UP and DOWN spin parts more thoroughly
% 
%---------------------------note 01.05.2013--------------------------------
% one reason for the non hermiticity of the hubbard hamilton matrix could 
% be the 2x2 lattice size. because than in the hopping process same states
% are created twice. sometimes even with different phases! run a hermitian
% check for 4x4 lattice and look if it still happens there!
%
% -----------------------IMPORTAN NOTE 08.05.13----------------------------
% rewrite this shit for 1D(and maybe 3D) and non half filling and with 
% option to exlude point group symmetries
% -------------------------------------------------------------------------
%------------------NOTE: 22.05.13 (ad basisCombination)--------------------
% METHOD: similar to Fehske, Weisse: (gives me wrong number for 2x2 but 
% right for 3x3, and very expensive
% METHOD LaFlorience: gives me right number for 2x2 and 3x3 and is faster
% TODO check for differences in basis states and why fehske weisse method
% fails for 2x2 lattice( apperently weisse fehske just wrong for 2x2
%--------------------------------------------------------------------------
% NOTE: 27.05.13: run for 4x4 including PG syms took ~6,5h and gave
% incorrect number of representatives. there was no hamiltonian calculation
% but basis alone gave no memory problems yet.
% 
% NOTE: 11.06.13: YEAAAAAAAAAAAAAAAAAAAAAAH it wOrKsyzz! 
% changes for it to work:
% * changed inclusion of symmetry eigenvalues px etc. in calculation of 
% state phases in cubicSymmetry function! because this way i applied it
% doubly in hubbard state norm calculation for example and in every further
% application of cubicSymmetry().  but should only be applied once ( except
% in calcSymOpsMatrixOverlap() where i need it from symOp2Repr part)
% * point group symmetries depending on k-values correctly chosen now
% according to sandvik paper. although diagonal mirroring for k= (kx,kx)
% different from all the others! ( k vector lies on mirror plane for Px and
% Py eg. but not for Pd) 
% * implementation of LaFlorience basis reduction now correct with right
% symmetries! could calculate the norm of hubbard states in laFlorience
% part for speed reasons(TODO) 
% * right norm calculation without doubly symmetry eigenvalue inclusion!
% cleaned up phases inclusion in hamiltonian calculations to correctly
% include symmetry eigenvalues 
% * correct calculation of symmetry matrix overlap calculation: i do it
% now: apply inverse(-signs for translations and inversed order) on basis
% state two : <D| <-inverse ; so have to first apply inverse of sym ops and
% then also symmetry operation which brings UP state to representative! 
% AND here i have to include the symmetry eigenvalue of the operation
% bringing the UP spin to its representative!! 
% thats all folks!
% NOTE: 14.06.13: for nSpinUp = 14, nSpinDown = 5 warning:-----------------
% Warning: Integer operands are required for colon operator when used as index 
% > In HubbardPhase>@(x,y)sum(replBasisStates(x,indUpToSiteIndex(1:y-1,x)),2) at 72
%   In HubbardPhase at 72
%   In calcHamiltonDownPointGroup at 114
%   In calcHubbardHamiltonian at 115
%   In Hubbard_ED at 269 
% !!!---------------------------------------------------------------------
% NOTE: 24.06.13: begin to change to implement other lattices..
% start with rectangular
% NOTE: 10.07.13: changed implementation of lattice size for rectanguklar
% lattices [nSitesY;nSitesX] to make matrix like matlab row column type
% NOTE: 10.07.13: seems right now for rectangular cases too. still not sure
% why Pd symmetry is different to Px and Py application?!TODO
% NOTE: 13.07.13:  include other lattice types: start with 'quadratic' 
% triangular lattice only translation!

%% preperation
clear all; close all; rehash; renameMatlabWindow;
% get computer info: only for unix for now:
if isunix
    [~,cpuInfos] = unix('cat /proc/cpuinfo');
    moreInfos = evalc('configinfo');
    baseSaveDir = '/afs/itp.tugraz.at/user/dobrautz/Dropbox/Master/tmp_diplArbeit/Diplomarbeit/ResultsHubbardED/';
else 
    cpuInfos =[]; 
    moreInfos = [];
    baseSaveDir = 'D:\Program Files (x86)\Dropbox\Dropbox\Master\tmp_diplArbeit\Diplomarbeit\ResultsHubbardED\';
end

%% Inputparameters
paramU = 2;                         % on-site-repulsion parameter U
paramT = 1;                         % hopping parameter t
latticeSize = [10];                % lattice parameters
diffQuarterFillingUp = 0;        % diff. of spins to quarter filling UP
diffQuarterFillingDown = -0;        % diff. of spins to quarter filling DOWN
nLanczosEigenvalues = 1;            % lanczos eigenvalues
lanzcosPrecision = 1e-8;            % lanczos genauigkeit
lanczosSteps = 200;                 % max. lanczos schritte
flagApplyPointGroup = 0;            % flag for point group symmetry
flagDoNoSymmetry = 1;               % flag for no symmetry calculation
flagDoOnlyTrans = 0;                % flag for only translation symmetry
flagSaveResults = 0;                % flag for saving
flagDispTables = 1;              % flag for displaying table of results
latticeType = 'cubic';              % string indicating lattice Type
latticeSaveDir = ['1D',latticeType];% standard lattice Dim 1D! changed below
disp(['!!!!!!!!!!!     calculating: ',num2str(diffQuarterFillingUp),', ',num2str(diffQuarterFillingDown),'      !!!!!!!!'])
t_transcorr_on_site = 0; 
transcorr_param = 0.0; %[linspace(-2,2,100)]; 

%% Parameters
[indNeighbors, nSites] = getNearestNeighbors(latticeType,latticeSize);

[latticeDim,~] = size(latticeSize); % Dimension of lattice
oddSiteFactor = mod(nSites,2)/2;    % odd number of sites factor

% fillings for UP and DOWN spins dep. if even or odd lattice Sites
fillingUp = diffQuarterFillingUp + oddSiteFactor + nSites/2;
fillingDown = diffQuarterFillingDown + oddSiteFactor + nSites/2;
% number of basis states for different spin channel
nStatesUp = nchoosek(nSites,fillingUp);     
nStatesDown = nchoosek(nSites,fillingDown); 

% linear lattice length for k vector calculation, since always cubic:
linearLatticeSizeY = latticeSize(1);
% k-Vectors
kVector = 2*pi/linearLatticeSizeY * [0:linearLatticeSizeY-1];%#ok
% kVector = 2*pi/linearLatticeSize * [-linearLatticeSize/2:linearLatticeSize/2-1];%#ok
if latticeDim == 2
    linearLatticeSizeX = latticeSize(2);
    kVectorX = 2*pi/linearLatticeSizeX * [0:linearLatticeSizeX-1];%#ok

    [ky1,kx1] = ndgrid(kVector,kVectorX); kx1 = kx1(:)'; ky1 = ky1(:)';
    kVector = [kx1;ky1];
    
    % if 2D change latticeSave String
    latticeSaveDir = [num2str(latticeSize(1)),'X',num2str(latticeSize(2)),...
        latticeType];
end

disp(['N_u = ',num2str(fillingUp)])
disp(['N_d = ',num2str(fillingDown)])
disp(['size Hilbert space: ',num2str(nStatesUp*nStatesDown)])

disp(['memory needed for 2 Lanczos vectors: ',Bytes2str(2*8*nStatesUp*nStatesDown)])
disp(['worst case sparse Hamilton memory: ',Bytes2str(8*2*nStatesUp*nStatesDown*fillingUp*fillingDown*2)])
% keyboard
%% Basis Creation
% without translational and point group symmetries
% NOTE: need integer values of whole UP spin basis for right linking in UP 
% spin hopping part calculation of hamiltonian
[basisStatesUp,intStatesUp] = createHeisenbergBasis(nStatesUp, nSites, fillingUp);
basisStatesDown = createHeisenbergBasis(nStatesDown, nSites, fillingDown);

if flagDoOnlyTrans == 1% && nSites < 16
%% ONLY TRANSLATIONAL SYMMETRY (1D and 2D)
tStartTrans = tic;
% calculate representatives, symmetry operations leaving repr. invariant 
% and linking list and corresponding symmetry operations to repr.:
[basisReprUp, symOpInvariantsUp, index2ReprUp, symOp2ReprUp ] = ...
    findReprOnlyTrans( basisStatesUp, latticeSize );

% initialize container for results and count of states t:
resultsTransOnly = {'lattice size','U/t','# spin up', '# spin down','time','PC info';...
    latticeSize,paramU/paramT,sum(basisStatesUp(1,:),2),sum(basisStatesDown(1,:),2),[],[cpuInfos,moreInfos];
    'kx','ky','E0','# states',[],[]};

% loop over k-vector
for ikVector = 1:length(kVector)
    
    kValue = kVector(:,ikVector);
    
    % combine up and down spin states to hubbard basis laFlorience style
    % calculate the norm of hubbard basis states: here also incompatible 
    % representatives per kValue are sorted out
    [compDownStatesPerRepr, compInd2ReprDown,normHubbardStates] = combine2HubbardBasisOnlyTrans(...
        symOpInvariantsUp, basisStatesDown, latticeSize, kValue);
    
    % calculate hamiltonian:
    H = calcHubbardHamiltonian(basisReprUp, compDownStatesPerRepr, ...
        compInd2ReprDown, paramU , paramT , indNeighbors, normHubbardStates,...
        symOpInvariantsUp, kValue, index2ReprUp, symOp2ReprUp, intStatesUp,...
        basisStatesDown,latticeSize);
    
    % lanczos
	E0 = lanczos_myversion(H,nLanczosEigenvalues,lanzcosPrecision,lanczosSteps);
    
    % this can only be done for small systems! 
    if t_transcorr_on_site
        
        % also do the un-correlated: 
        [orig_V, orig_E0] = eig(full(H));
        [orig_E0, min_ind] = min(diag(orig_E0)); 
        orig_V = orig_V(:,min_ind); 
        
        hf_coeff_orig = max(abs(orig_V));
        
        ind = 1;
        gs_vec_trans = zeros(size(H,1),length(transcorr_param));
        hf_coeff = zeros(length(transcorr_param));
        
        for loop_param = transcorr_param
       
            t_mat = get_trans_mat_on_site(H, paramU, loop_param ); 

            trans_H = expm(-t_mat) * H * expm(t_mat); 

            [trans_V, trans_E0] = eig(trans_H); 
            
            % store the groundstate vector: 
            [~,min_ind] = min(real(diag(trans_E0)));
            gs_vec_trans(:,ind) = sort(abs(trans_V(:,min_ind)),'descend'); 
            hf_coeff(ind) = gs_vec_trans(1,ind); 
            
            ind = ind + 1;
            
            if abs(orig_E0 - trans_E0(min_ind,min_ind)) > 1.0e-8
                
                print('eigenvalues differ!')
                
            end
            
        end
        
%         figure(ikVector) 
%         plot(sort(abs(orig_V),'descend'),'r*')
%         hold on
%         plot(gs_vec_trans)
%         hold off
        
        figure(ikVector)
        plot(0, hf_coeff_orig,'r*')
        hold on
        plot(transcorr_param, hf_coeff)
        xlabel('J')
        ylabel('HF Coeff.')
        title(['U = ', num2str(paramU), ', k = ', num2str(kValue)])
        hold off
        
    end
    % number of basis states per kValue
    specNumStates = sum(cellfun(@length,normHubbardStates));
    
    % save results
    resultsTransOnly{ikVector+3,3} = E0;
    resultsTransOnly{ikVector+3,4} = specNumStates;
    resultsTransOnly{ikVector+3,5} = [];
    resultsTransOnly{ikVector+3,6} = [];
    
    resultsTransOnly{ikVector+3,1} = kValue(1);
    resultsTransOnly{ikVector+3,2} = [];
    
    if latticeDim == 2
        resultsTransOnly{ikVector+3,2} = kValue(2);
    end
    
    clear H normHubbardStates comp* 

end

% save results:
timeOnlyTrans = toc(tStartTrans);

resultsTransOnly{2,5} = timeOnlyTrans; 

if flagSaveResults == 1
%     savedir = [baseSaveDir,'OnlyTranslation',filesep,latticeSaveDir,filesep];
    
    if latticeDim == 1
        savename = [date,'_',num2str(latticeSize(1)),'x','1','_U=',num2str(paramU),...
            '_t=',num2str(paramT),'_nUp=',num2str(sum(basisStatesUp(1,:),2)),...
            '_nDown=',num2str(sum(basisStatesDown(1,:),2)),'TransOnly'];
        
    else
        savename = [date,'_',num2str(latticeSize(1)),'x',num2str(latticeSize(2)),...
            '_U=',num2str(paramU), '_t=',num2str(paramT),'_nUp=',...
            num2str(sum(basisStatesUp(1,:),2)),'_nDown=',...
            num2str(sum(basisStatesDown(1,:),2)),'TransOnly'];
        
    end
    
%     save([savedir,savename,'.mat'],'resultsTransOnly');
    save([savename,'.mat'],'resultsTransOnly');
end

if flagDispTables == 1
    resultsTransOnly%#ok
    sum(cellfun(@sum,resultsTransOnly(cellfun(@isnumeric,resultsTransOnly(:,4)),4)))

end

end

if flagApplyPointGroup == 1 && latticeDim == 2
%% Two Dimensions point group and translations
% ---------------------- NOTE: 24.05.13: ----------------------------------
% now the change in the old code begins...
% difference now with point group symmetries: it depends on the k-vector
% value which lattice symmetry can be applied: so the calculation of the
% representatives also depends on the k vector
% see: Comp. Studies of Quant. Spin Systems, A.W. Sandvik p.117
% -------------------------------------------------------------------------

% initialize container for results 
resultsPointGroup = {'lattice size', 'U/t', '# spin UP', '# spin DOWN','time','PC info',[],[];...
    latticeSize,paramU/paramT,sum(basisStatesUp(1,:),2),sum(basisStatesDown(1,:),2),[],[cpuInfos,moreInfos],[],[];...
    'kx','ky','px','py','pd','pe','E0','# states'};

[~,pcName] = system('hostname');
handleWaitbar = waitbar(0,pcName);
set(handleWaitbar,'name',pcName);

tStartPoint = tic;
% loop over k-values (trans. symm.)
for ikVector = 1:numel(kx1)          
    
    kx = kVector(1,ikVector); 
    ky = kVector(2,ikVector);
    kValue = kVector(:,ikVector);
    %% Applicable lattice symmetries
    % see: Comp. Studies of Quant. Spin Systems, A.W. Sandvik p.117
    if (kx ~= 0 && ky == 0)%|| (abs(kx) == pi && abs(ky) ~= pi)
        for py = [-1,1]
            % Q = (1 + py*Py);
            % calculate representatives, symmetry operations leaving repr. invariant 
            % and linking list and corresponding symmetry operations to repr.:
            [basisReprUp, symOpInvariantsUp, index2ReprUp, symOp2ReprUp ] = ...
                PointGroupFindRepr( basisStatesUp, 'py', kValue,latticeSize, py );
            % combine up and down spin states to hubbard basis: defined by combination
            % of basisReprUp and downStatesPerRepr. Method LaFlorience:
            % also includes norm calculation and compatibility test
            [compDownStatesPerRepr,compInd2ReprDown,normHubbardStates] = ...
                PointGroupLaFlorience( symOpInvariantsUp, basisStatesDown, latticeSize);

            % calculate hamiltonian:
            H = calcHubbardHamiltonian(basisReprUp, compDownStatesPerRepr,...
                compInd2ReprDown,paramU, paramT, indNeighbors, normHubbardStates, ...
                symOpInvariantsUp, kValue, index2ReprUp, symOp2ReprUp, intStatesUp, ...
                basisStatesDown , latticeSize);
            
            % count sates
            specNumStates = sum(cellfun(@length,normHubbardStates));
            % lanczos
            E0 = lanczos_myversion(H,nLanczosEigenvalues,lanzcosPrecision,lanczosSteps);
            % table update
            resultsPointGroup = tableUpdate(resultsPointGroup,E0,specNumStates,kx,ky,py);
            
        end
        
    elseif (kx == 0 && ky ~= 0)  %|| (abs(kx) ~= pi && abs(ky) == pi)
        for px = [-1,1]
            % Q = (1 + px*Px);
            [basisReprUp, symOpInvariantsUp, index2ReprUp, symOp2ReprUp ] = ...
                PointGroupFindRepr( basisStatesUp, 'px',kValue,latticeSize, px );
            
            [compDownStatesPerRepr,compInd2ReprDown,normHubbardStates] = ...
                PointGroupLaFlorience( symOpInvariantsUp, basisStatesDown, latticeSize);
            
            H = calcHubbardHamiltonian(basisReprUp, compDownStatesPerRepr,...
                compInd2ReprDown,paramU, paramT, indNeighbors, normHubbardStates, ...
                symOpInvariantsUp, kValue, index2ReprUp, symOp2ReprUp, intStatesUp, ...
                basisStatesDown, latticeSize );

            specNumStates = sum(cellfun(@length,normHubbardStates));
            
            E0 = lanczos_myversion(H,nLanczosEigenvalues,lanzcosPrecision,lanczosSteps);
            
            resultsPointGroup = tableUpdate(resultsPointGroup,E0,specNumStates,kx,ky,px);
            
        end
        
    elseif (kx == ky && (kx ~= 0 && abs(kx) ~= pi)) && latticeSize(1)==latticeSize(2)
        for pd = [-1,1]
            % Q = (1 + pd*Pd);
            [basisReprUp, symOpInvariantsUp, index2ReprUp, symOp2ReprUp ] = ...
                PointGroupFindRepr( basisStatesUp, 'pd',kValue,latticeSize, pd );
            
            [compDownStatesPerRepr,compInd2ReprDown,normHubbardStates] = ...
                PointGroupLaFlorience( symOpInvariantsUp, basisStatesDown, latticeSize);

            H = calcHubbardHamiltonian(basisReprUp, compDownStatesPerRepr,...
                compInd2ReprDown,paramU, paramT, indNeighbors, normHubbardStates, ...
                symOpInvariantsUp, kValue, index2ReprUp, symOp2ReprUp, intStatesUp, ...
                basisStatesDown, latticeSize );

            specNumStates = sum(cellfun(@length,normHubbardStates));
            
            E0 = lanczos_myversion(H,nLanczosEigenvalues,lanzcosPrecision,lanczosSteps);
            
            resultsPointGroup = tableUpdate(resultsPointGroup,E0,specNumStates,kx,ky,pd);
        end
        
    elseif (kx == -ky && kx ~= 0) && latticeSize(1)==latticeSize(2)
        for pe = [-1,1]
            % Q = (1 + pe*Pe);
            [basisReprUp, symOpInvariantsUp, index2ReprUp, symOp2ReprUp ] = ...
                PointGroupFindRepr( basisStatesUp, 'pe',kValue,latticeSize, pe );
            
            [compDownStatesPerRepr,compInd2ReprDown,normHubbardStates] = ...
                PointGroupLaFlorience( symOpInvariantsUp, basisStatesDown, latticeSize);
            
            H = calcHubbardHamiltonian(basisReprUp, compDownStatesPerRepr,...
                compInd2ReprDown,paramU, paramT, indNeighbors, normHubbardStates, ...
                symOpInvariantsUp, kValue, index2ReprUp, symOp2ReprUp, intStatesUp, ...
                basisStatesDown, latticeSize );
            
            specNumStates = sum(cellfun(@length,normHubbardStates));
            
            E0 = lanczos_myversion(H,nLanczosEigenvalues,lanzcosPrecision,lanczosSteps);
            
            resultsPointGroup = tableUpdate(resultsPointGroup,E0,specNumStates,kx,ky,pe);
        end
        
    elseif (kx == ky && (kx == 0 || abs(kx) == pi))
        for px = [-1,1]
            for py = [-1,1]
                if px == py && latticeSize(1) == latticeSize(2)
                    for pd = [-1,1]
                        % Q = (1 + pd*Pd)*(1 + py*Py)*(1 + px*Px);
                        [basisReprUp, symOpInvariantsUp, index2ReprUp, symOp2ReprUp ] = ...
                            PointGroupFindRepr( basisStatesUp, '',kValue,latticeSize, pd, py, px );
                        
                        [compDownStatesPerRepr,compInd2ReprDown,normHubbardStates] = ...
                            PointGroupLaFlorience( symOpInvariantsUp, basisStatesDown, latticeSize);
                        
                        H = calcHubbardHamiltonian(basisReprUp, compDownStatesPerRepr,...
                            compInd2ReprDown,paramU, paramT, indNeighbors, normHubbardStates, ...
                            symOpInvariantsUp, kValue, index2ReprUp, symOp2ReprUp, intStatesUp, ...
                            basisStatesDown, latticeSize );

                        specNumStates = sum(cellfun(@length,normHubbardStates));
                        
                        [E0,V0] = lanczos_myversion(H,nLanczosEigenvalues,lanzcosPrecision,lanczosSteps);
                        
                        resultsPointGroup = tableUpdate(resultsPointGroup,E0,specNumStates,kx,ky,px,py,pd);
                    end
                else
                    % Q = (1 + py*Py)*(1 + px*Px);
                    [basisReprUp, symOpInvariantsUp, index2ReprUp, symOp2ReprUp ] = ...
                            PointGroupFindRepr( basisStatesUp, '',kValue,latticeSize, py, px );
                    
                    [compDownStatesPerRepr,compInd2ReprDown,normHubbardStates] = ...
                        PointGroupLaFlorience( symOpInvariantsUp, basisStatesDown, latticeSize); 

                    H = calcHubbardHamiltonian(basisReprUp, compDownStatesPerRepr,...
                        compInd2ReprDown,paramU, paramT, indNeighbors, normHubbardStates, ...
                        symOpInvariantsUp, kValue, index2ReprUp, symOp2ReprUp, intStatesUp, ...
                        basisStatesDown, latticeSize );
  
                    specNumStates = sum(cellfun(@length,normHubbardStates));
                    
                    E0 = lanczos_myversion(H,nLanczosEigenvalues,lanzcosPrecision,lanczosSteps);
                    
                    resultsPointGroup = tableUpdate(resultsPointGroup,E0,specNumStates,kx,ky,px,py);
                end
            end
        end
        
    else
        % Q = 1;
        [basisReprUp, symOpInvariantsUp, index2ReprUp, symOp2ReprUp ] = ...
            PointGroupFindRepr( basisStatesUp, '',kValue,latticeSize);
        
       [compDownStatesPerRepr,compInd2ReprDown,normHubbardStates] = ...
            PointGroupLaFlorience( symOpInvariantsUp, basisStatesDown, latticeSize);
        
        H = calcHubbardHamiltonian(basisReprUp, compDownStatesPerRepr,...
            compInd2ReprDown,paramU, paramT, indNeighbors, normHubbardStates, ...
            symOpInvariantsUp, kValue, index2ReprUp, symOp2ReprUp, intStatesUp, ...
            basisStatesDown, latticeSize );
        
        specNumStates = sum(cellfun(@length,normHubbardStates));
        
        E0 = lanczos_myversion(H,nLanczosEigenvalues,lanzcosPrecision,lanczosSteps);
        
        resultsPointGroup = tableUpdate(resultsPointGroup,E0,specNumStates,kx,ky);
    end
    
    clear H comp* norm*
    
    waitbar(ikVector/numel(kx1),handleWaitbar,['k = ',num2str(ikVector),' of ',num2str(numel(kx1))])
end

timePointGroup = toc(tStartPoint);

close(handleWaitbar)


% save results
resultsPointGroup{2,5} = timePointGroup;

if flagSaveResults == 1
%     savedir = [baseSaveDir,'PointGroup',filesep,latticeSaveDir,filesep];
    
    savename = [date,'_',num2str(latticeSize(1)),'x',num2str(latticeSize(2)),...
        '_U=',num2str(paramU),'_t=',num2str(paramT),'_nUp=',num2str(sum(basisStatesUp(1,:),2)),...
        '_nDown=',num2str(sum(basisStatesDown(1,:),2)),'PointGroup'];
    
    
%     save([savedir,savename,'.mat'],'resultsPointGroup');
    save([savename,'.mat'],'resultsPointGroup');
end

if flagDispTables == 1
    resultsPointGroup%#ok
    
sum(cellfun(@sum,resultsPointGroup(cellfun(@isnumeric,resultsPointGroup(:,end)),end)))

end

end



if flagDoNoSymmetry == 1% && nSites < 10
%% NO Symmetry case:
    tStartNoSym = tic;
    
    % calc hamiltonian
    HamNosym = HubbardHamiltonNoSymmetries(basisStatesUp, basisStatesDown, ...
        indNeighbors, paramU, paramT);
    
    % lanzcos
    E0_full = lanczos_myversion(HamNosym,nLanczosEigenvalues,lanzcosPrecision,lanczosSteps);
    
    timeNoSymmetry = toc(tStartNoSym);
    
     % this can only be done for small systems! 
    if t_transcorr_on_site
        
        % also do the un-correlated: 
        [orig_V, orig_E0] = eig(full(HamNosym));
        [orig_E0, min_ind] = min(real(diag(orig_E0)));
        orig_V = orig_V(:,min_ind);
        hf_coeff_orig = max(abs(orig_V));
        
        ind = 1;
        gs_vec_trans = zeros(size(HamNosym,1),length(transcorr_param));
        gs_vec_trans_next = zeros(size(HamNosym,1),length(transcorr_param));
        gs_vec_trans_hop = zeros(size(HamNosym,1), length(transcorr_param));
        
        % also do some combinations.. but i need to decide 
        % on the sign and on the relations between the 
        % parameters.. no.. since not all of the transcorr factors 
        % commute i cannot just simply add them.. i mean here i can, 
        % since the systems are so small, but in principle I cant
        
        gs_vec_trans_sum = zeros(size(HamNosym,1), length(transcorr_param));
        
        t_mat = get_trans_mat_on_site(HamNosym, paramU, 1.0 ); 
        t_mat_next = get_trans_mat_neighbor(basisStatesUp, basisStatesDown, ... 
                indNeighbors, 1.0); 
        t_mat_hop = get_trans_mat_hopping(HamNosym, 1.0);
        t_mat_sum = t_mat + t_mat_next; 
        
        hf_coeff_next = zeros(length(transcorr_param),1); 
        hf_coeff = zeros(length(transcorr_param),1); 
        hf_coeff_hop = zeros(length(transcorr_param),1);
        hf_coeff_sum = zeros(length(transcorr_param),1);
        
        for loop_param = transcorr_param
            
            trans_H = expm(-loop_param * t_mat) * HamNosym * expm(loop_param * t_mat); 
            
            trans_H_next = expm(-loop_param * t_mat_next) * HamNosym * expm(loop_param * t_mat_next); 
            
            trans_H_hop = expm(-loop_param * t_mat_hop) * HamNosym * expm(loop_param * t_mat_hop);
            
            trans_H_sum = expm(-loop_param * t_mat_sum) * HamNosym * expm(loop_param * t_mat_sum);
            [trans_V, trans_E0] = eig(trans_H); 
            [trans_V_next, trans_E0_next] = eig(trans_H_next); 
            [trans_V_hop, trans_E0_hop] = eig(trans_H_hop);
            [trans_V_sum, trans_E0_sum] = eig(trans_H_sum);
            
            % store the groundstate vector: 
            [trans_E0_next,min_ind] = min(real(diag(trans_E0_next)));
            gs_vec_trans_next(:,ind) = sort(abs(trans_V_next(:,min_ind)),'descend'); 
            hf_coeff_next(ind) = gs_vec_trans_next(1,ind); 

            [trans_E0,min_ind] = min(real(diag(trans_E0)));
            gs_vec_trans(:,ind) = sort(abs(trans_V(:,min_ind)),'descend'); 
            hf_coeff(ind) = gs_vec_trans(1,ind); 

            [trans_E0_hop, min_ind] = min(real(diag(trans_E0_hop)));
            gs_vec_trans_hop(:,ind) = sort(abs(trans_V_hop(:,min_ind)), 'descend');
            hf_coeff_hop(ind) = gs_vec_trans_hop(1,ind); 
            
            [trans_E0_sum, min_ind] = min(real(diag(trans_E0_sum))); 
            gs_vec_trans_sum(:,ind) = sort(abs(trans_V_sum(:,min_ind)), 'descend');
            hf_coeff_sum(ind) = gs_vec_trans_sum(1,ind); 
            
            disp(['iteration ', num2str(ind), ' done'])
            ind = ind + 1;
            
            if abs(orig_E0 - trans_E0) > 1.0e-8
                
                disp(' on-site eigenvalues differ!')
                
            end
            if abs(orig_E0 - trans_E0_next) > 1.0e-8
                
                disp('neighbor eigenvalues differ!')
                
            end
            if abs(orig_E0 - trans_E0_hop) > 1.0e-8
                disp('hopping eigenvalues differ!')
            end 
            if abs(orig_E0 - trans_E0_sum) > 1.0e-8
                disp('sum eigenvalues differ!')
            end
        end
        
%         figure(1) 
%         plot(sort(abs(orig_V),'descend'),'r*')
%         hold on
%         plot(gs_vec_trans)
%         hold off
        
        figure(length(kVector)+1)
        plot(0, hf_coeff_orig,'r*')
        hold on
        plot(transcorr_param, hf_coeff, transcorr_param, hf_coeff_next, ...
            transcorr_param, hf_coeff_hop, transcorr_param, hf_coeff_sum)
        legend('orig. Hamiltonian','on-site corr.', 'neighbor corr.','hop corr.', 'sum')
        xlabel('J')
        ylabel('HF Coeff.')
        title(['nOrbs: ', num2str(latticeSize),', U: ', num2str(paramU), ...
            ', nel: ',num2str(fillingDown+fillingUp)])
        
    end
    
    % save results
    resultsNoSymmetry = {'lattice size', 'U/t', '# spin UP', '# spin DOWN',...
        'time','E0','# states','PC info';  latticeSize, paramU/paramT, sum(basisStatesUp(1,:),2),...
        sum(basisStatesDown(1,:),2), timeNoSymmetry,E0_full,nStatesUp*nStatesDown,[cpuInfos,moreInfos]};


    if flagSaveResults == 1
        savedir = [baseSaveDir,'noSymmetry',filesep,latticeSaveDir,filesep];
        
        if latticeDim == 1
            savename = [date,'_',num2str(latticeSize(1)),'x','1','_U=',num2str(paramU),...
                '_t=',num2str(paramT),'_nUp=',num2str(sum(basisStatesUp(1,:),2)),...
                '_nDown=',num2str(sum(basisStatesDown(1,:),2)),'noSymmetry'];
        else
            savename = [date,'_',num2str(latticeSize(1)),'x',num2str(latticeSize(2)),'_U=',num2str(paramU),...
                '_t=',num2str(paramT),'_nUp=',num2str(sum(basisStatesUp(1,:),2)),...
                '_nDown=',num2str(sum(basisStatesDown(1,:),2)),'noSymmetry'];
        end
        
        save([savedir,savename,'.mat'],'resultsNoSymmetry');
    end
    
    if flagDispTables == 1;
        resultsNoSymmetry%#ok
    end
    
    
end

if flagSaveResults == 1
    clear all
    exit
end
