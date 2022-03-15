function callHubbard_ED
%callHubbard_ED function version of Hubbard_ED script for compiling

%% preperation
% check if matlab pool is open
if matlabpool('size') > 0
    matlabpool close;
end
% open matlabpool with 'local' setting, check for licenses
matlabpool;

% get computer info: only for unix for now:
if isunix
    [~,cpuInfos] = unix('cat /proc/cpuinfo');
    moreInfos = evalc('configinfo');
else
    cpuInfos =[];
    moreInfos = [];
end

%% Inputparameters
% load and assign input data
if isunix
    cd('/afs/itp.tugraz.at/user/dobrautz/Dropbox/Master/tmp_diplArbeit/Diplomarbeit/tags/stableHubbardED_eigenValueCalc17Juni2013');
    inputData = load('/afs/itp.tugraz.at/user/dobrautz/Dropbox/Master/tmp_diplArbeit/Diplomarbeit/tags/stableHubbardED_eigenValueCalc17Juni2013/inputHubbard_ED.txt');
    system('export MCR_CACHE_ROOT=/tmp/dobrautz')
else
    cd('D:\Program Files (x86)\Dropbox\Dropbox\Master\tmp_diplArbeit\Diplomarbeit\tags\stableHubbardED_eigenValueCalc17Juni2013');
    inputData = load('D:\Program Files (x86)\Dropbox\Dropbox\Master\tmp_diplArbeit\Diplomarbeit\tags\stableHubbardED_eigenValueCalc17Juni2013\inputHubbard_ED.txt');
end

paramU = inputData(1);                     % on-site-repulsion parameter U
paramT = inputData(2);                     % hopping parameter t
latticeSize = [inputData(3);inputData(4)]; % lattice parameters
diffQuarterFillingUp = inputData(5);       % diff. of spins to quarter filling UP
diffQuarterFillingDown = inputData(6);     % diff. of spins to quarter filling DOWN
nLanczosEigenvalues = inputData(7);        % lanczos eigenvalues
lanczosPrecision = inputData(8);           % lanczos genauigkeit
lanczosSteps = inputData(9);               % max. lanczos schritte
flagApplyPointGroup = inputData(10);       % flag for point group symmetry
flagDoNoSymmetry = inputData(11);          % flag for no symmetry calculation
flagDoOnlyTrans = inputData(12);           % flag for only translation symmetry
flagSaveResults = inputData(13);           % flag for saving

%% Parameters
nSites = prod(latticeSize);         % number of sites
[latticeDim,~] = size(latticeSize); % Dimension of lattice
oddSiteFactor = mod(nSites,2)/2;    % odd number of sites factor
% fillings for UP and DOWN spins dep. if even or odd lattice Sites
fillingUp = diffQuarterFillingUp + oddSiteFactor;
fillingDown = diffQuarterFillingDown + oddSiteFactor;
% number of basis states for different spin channel
nStatesUp = nchoosek(nSites,fillingUp + nSites/2);      
nStatesDown = nchoosek(nSites,fillingDown + nSites/2); 
% Indices of nearest neighbors (is dimension indication in later functions)
indNeighbors = findCubicNearestNeighbors(latticeSize);
% linear lattice length for k vector calculation, since always cubic:
linearLatticeSize = latticeSize(1);   
% k-Vectors
kVector = 2*pi/linearLatticeSize * [0:linearLatticeSize-1];%#ok
% kVector = 2*pi/linearLatticeSize * [-linearLatticeSize/2:linearLatticeSize/2-1];%#ok
if latticeDim == 2
    [kx1,ky1] = ndgrid(kVector); kx1 = kx1(:)'; ky1 = ky1(:)';
    kVector = [kx1;ky1];
end

%% Basis Creation
% without translational and point group symmetries
% NOTE: need integer values of whole UP spin basis for right linking in UP 
% spin hopping part calculation of hamiltonian
[basisStatesUp,intStatesUp] = createHeisenbergBasis(nStatesUp, nSites, fillingUp);
basisStatesDown = createHeisenbergBasis(nStatesDown, nSites, fillingDown);

if flagDoOnlyTrans == 1% && nSites < 16
%% ONLY TRANSLATIONAL SYMMETRY (1D and 2D)
tic
% calculate representatives, symmetry operations leaving repr. invariant 
% and linking list and corresponding symmetry operations to repr.:
[basisReprUp, symOpInvariantsUp, index2ReprUp, symOp2ReprUp ] = ...
    findReprOnlyTrans( basisStatesUp, latticeDim );

% inittialize container for results and count of states t:
resultsTransOnly = {'lattice size','U/t','# spin up', '# spin down','time','PC info';...
    latticeSize,paramU/paramT,sum(basisStatesUp(1,:),2),sum(basisStatesDown(1,:),2),[],[cpuInfos,moreInfos];
    'kx','ky','E0','# states',[],[]};


% loop over k-vector
for ikVector = 1:length(kVector)
    
    kValue = kVector(:,ikVector);
    
    % calculate the norm of hubbard basis states: here also incompatible 
    % representatives per kValue are sorted out
    % combine up and down spin states to hubbard basis: defined by combination
    % of basisReprUp and downStatesPerRepr
    [compDownStatesPerRepr, compInd2ReprDown,normHubbardStates] = combine2HubbardBasisOnlyTrans(...
        symOpInvariantsUp, basisStatesDown, latticeDim, kValue);
    
    % calculate hamiltonian:
    H = calcHubbardHamiltonian(basisReprUp, compDownStatesPerRepr, ...
        compInd2ReprDown, paramU , paramT , indNeighbors, normHubbardStates,...
        symOpInvariantsUp, kValue, index2ReprUp, symOp2ReprUp, intStatesUp,...
        basisStatesDown);
    
    % lanczos
	E0 = lanczos_myversion(H,nLanczosEigenvalues,lanczosPrecision,lanczosSteps);
    
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

timeOnlyTrans = toc;

resultsTransOnly{2,5} = timeOnlyTrans; %#ok

if flagSaveResults == 1 
    
    if isunix
        savedir = ['/afs/itp.tugraz.at/user/dobrautz/Dropbox/Master/Diplomarbeit/meineProgramme/Hubbard_ED/Results/OnlyTranslation',filesep];
    else
        savedir = ['D:\Program Files (x86)\Dropbox\Dropbox\Master\Diplomarbeit\meineProgramme\Hubbard_ED\Results\OnlyTranslation',filesep];
    end
    
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
    
    save([savedir,savename,'.mat'],'resultsTransOnly');
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
tic
% initialize container for results 
resultsPointGroup = {'lattice size', 'U/t', '# spin UP', '# spin DOWN','time','PC info',[],[];...
    latticeSize,paramU/paramT,sum(basisStatesUp(1,:),2),sum(basisStatesDown(1,:),2),[],[cpuInfos,moreInfos],[],[];...
    'kx','ky','px','py','pd','pe','E0','# states'};

% loop over k-values (trans. symm.)
for ikVector = 1:numel(kx1)          
    
    kx = kVector(1,ikVector); 
    ky = kVector(2,ikVector);
    kValue = kVector(:,ikVector);
    %% Applicable lattice symmetries
    % see: Comp. Studies of Quant. Spin Systems, A.W. Sandvik p.117
    if  (kx ~= 0 && ky == 0) %|| (abs(kx) == pi && abs(ky) ~= pi)
        for py = [-1,1]
            % Q = (1 + py*Py);
            % calculate representatives, symmetry operations leaving repr. invariant 
            % and linking list and corresponding symmetry operations to repr.:
            [basisReprUp, symOpInvariantsUp, index2ReprUp, symOp2ReprUp ] = ...
                PointGroupFindRepr( basisStatesUp, 'py', kValue, py );
            % combine up and down spin states to hubbard basis: defined by combination
            % of basisReprUp and downStatesPerRepr. Method LaFlorience:
            % also includes norm calculation and compatibility test
            [compDownStatesPerRepr,compInd2ReprDown,normHubbardStates] = ...
                PointGroupLaFlorience( symOpInvariantsUp, basisStatesDown);

            % calculate hamiltonian:
            H = calcHubbardHamiltonian(basisReprUp, compDownStatesPerRepr,...
                compInd2ReprDown,paramU, paramT, indNeighbors, normHubbardStates, ...
                symOpInvariantsUp, kValue, index2ReprUp, symOp2ReprUp, intStatesUp, ...
                basisStatesDown );
            
            % count sates
            specNumStates = sum(cellfun(@length,normHubbardStates));
            % lanczos
            E0 = lanczos_myversion(H,nLanczosEigenvalues,lanczosPrecision,lanczosSteps);
            % table update
            resultsPointGroup = tableUpdate(resultsPointGroup,E0,specNumStates,kx,ky,py);
            
        end
        
    elseif (kx == 0 && ky ~= 0) %|| (abs(kx) ~= pi && abs(ky) == pi)
        for px = [-1,1]
            % Q = (1 + px*Px);
            [basisReprUp, symOpInvariantsUp, index2ReprUp, symOp2ReprUp ] = ...
                PointGroupFindRepr( basisStatesUp, 'px',kValue, px );
            
            [compDownStatesPerRepr,compInd2ReprDown,normHubbardStates] = ...
                PointGroupLaFlorience( symOpInvariantsUp, basisStatesDown);
            
            H = calcHubbardHamiltonian(basisReprUp, compDownStatesPerRepr,...
                compInd2ReprDown,paramU, paramT, indNeighbors, normHubbardStates, ...
                symOpInvariantsUp, kValue, index2ReprUp, symOp2ReprUp, intStatesUp, ...
                basisStatesDown );

            specNumStates = sum(cellfun(@length,normHubbardStates));
            
            E0 = lanczos_myversion(H,nLanczosEigenvalues,lanczosPrecision,lanczosSteps);
            
            resultsPointGroup = tableUpdate(resultsPointGroup,E0,specNumStates,kx,ky,px);
            
        end
        
    elseif (kx == ky && (kx ~= 0 && abs(kx) ~= pi))
        for pd = [-1,1]
            % Q = (1 + pd*Pd);
            [basisReprUp, symOpInvariantsUp, index2ReprUp, symOp2ReprUp ] = ...
                PointGroupFindRepr( basisStatesUp, 'pd',kValue, pd );
            
            [compDownStatesPerRepr,compInd2ReprDown,normHubbardStates] = ...
                PointGroupLaFlorience( symOpInvariantsUp, basisStatesDown);

            H = calcHubbardHamiltonian(basisReprUp, compDownStatesPerRepr,...
                compInd2ReprDown,paramU, paramT, indNeighbors, normHubbardStates, ...
                symOpInvariantsUp, kValue, index2ReprUp, symOp2ReprUp, intStatesUp, ...
                basisStatesDown );

            specNumStates = sum(cellfun(@length,normHubbardStates));
            
            E0 = lanczos_myversion(H,nLanczosEigenvalues,lanczosPrecision,lanczosSteps);
            
            resultsPointGroup = tableUpdate(resultsPointGroup,E0,specNumStates,kx,ky,pd);
        end
        
    elseif (kx == -ky && kx ~= 0)
        for pe = [-1,1]
            % Q = (1 + pe*Pe);
            [basisReprUp, symOpInvariantsUp, index2ReprUp, symOp2ReprUp ] = ...
                PointGroupFindRepr( basisStatesUp, 'pe',kValue, pe );
            
            [compDownStatesPerRepr,compInd2ReprDown,normHubbardStates] = ...
                PointGroupLaFlorience( symOpInvariantsUp, basisStatesDown);
            
            H = calcHubbardHamiltonian(basisReprUp, compDownStatesPerRepr,...
                compInd2ReprDown,paramU, paramT, indNeighbors, normHubbardStates, ...
                symOpInvariantsUp, kValue, index2ReprUp, symOp2ReprUp, intStatesUp, ...
                basisStatesDown );
            
            specNumStates = sum(cellfun(@length,normHubbardStates));
            
            E0 = lanczos_myversion(H,nLanczosEigenvalues,lanczosPrecision,lanczosSteps);
            
            resultsPointGroup = tableUpdate(resultsPointGroup,E0,specNumStates,kx,ky,pe);
        end
        
    elseif (kx == ky && (kx == 0 || abs(kx) == pi))
        for px = [-1,1]
            for py = [-1,1]
                if px == py%  && px==2
                    for pd = [-1,1]
                        % Q = (1 + pd*Pd)*(1 + py*Py)*(1 + px*Px);

                        [basisReprUp, symOpInvariantsUp, index2ReprUp, symOp2ReprUp ] = ...
                            PointGroupFindRepr( basisStatesUp, '',kValue, pd, py, px );
                        
                        [compDownStatesPerRepr,compInd2ReprDown,normHubbardStates] = ...
                            PointGroupLaFlorience( symOpInvariantsUp, basisStatesDown);
                        
                        H = calcHubbardHamiltonian(basisReprUp, compDownStatesPerRepr,...
                            compInd2ReprDown,paramU, paramT, indNeighbors, normHubbardStates, ...
                            symOpInvariantsUp, kValue, index2ReprUp, symOp2ReprUp, intStatesUp, ...
                            basisStatesDown );

                        specNumStates = sum(cellfun(@length,normHubbardStates));
                        
                        E0 = lanczos_myversion(H,nLanczosEigenvalues,lanczosPrecision,lanczosSteps);
                        
                        resultsPointGroup = tableUpdate(resultsPointGroup,E0,specNumStates,kx,ky,px,py,pd);
                    end
                else
                    % Q = (1 + py*Py)*(1 + px*Px);
                    [basisReprUp, symOpInvariantsUp, index2ReprUp, symOp2ReprUp ] = ...
                            PointGroupFindRepr( basisStatesUp, '',kValue, py, px );
                    
                    [compDownStatesPerRepr,compInd2ReprDown,normHubbardStates] = ...
                        PointGroupLaFlorience( symOpInvariantsUp, basisStatesDown); 

                    H = calcHubbardHamiltonian(basisReprUp, compDownStatesPerRepr,...
                        compInd2ReprDown,paramU, paramT, indNeighbors, normHubbardStates, ...
                        symOpInvariantsUp, kValue, index2ReprUp, symOp2ReprUp, intStatesUp, ...
                        basisStatesDown );
  
                    specNumStates = sum(cellfun(@length,normHubbardStates));
                    
                    E0 = lanczos_myversion(H,nLanczosEigenvalues,lanczosPrecision,lanczosSteps);
                    
                    resultsPointGroup = tableUpdate(resultsPointGroup,E0,specNumStates,kx,ky,px,py);
                end
            end
        end
        
    else
        % Q = 1;
        [basisReprUp, symOpInvariantsUp, index2ReprUp, symOp2ReprUp ] = ...
            PointGroupFindRepr( basisStatesUp, '',kValue);
        
       [compDownStatesPerRepr,compInd2ReprDown,normHubbardStates] = ...
            PointGroupLaFlorience( symOpInvariantsUp, basisStatesDown);
        
        H = calcHubbardHamiltonian(basisReprUp, compDownStatesPerRepr,...
            compInd2ReprDown,paramU, paramT, indNeighbors, normHubbardStates, ...
            symOpInvariantsUp, kValue, index2ReprUp, symOp2ReprUp, intStatesUp, ...
            basisStatesDown );
        
        specNumStates = sum(cellfun(@length,normHubbardStates));
        
        E0 = lanczos_myversion(H,nLanczosEigenvalues,lanczosPrecision,lanczosSteps);
        
        resultsPointGroup = tableUpdate(resultsPointGroup,E0,specNumStates,kx,ky);
    end
    
    clear H comp* norm*
    
end

timePointGroup = toc;

resultsPointGroup{2,5} = timePointGroup;%#ok

if flagSaveResults == 1
    if isunix
        savedir = ['/afs/itp.tugraz.at/user/dobrautz/Dropbox/Master/Diplomarbeit/meineProgramme/Hubbard_ED/Results/PointGroup',filesep];
    else
        savedir = ['D:\Program Files (x86)\Dropbox\Dropbox\Master\Diplomarbeit\meineProgramme\Hubbard_ED\Results\PointGroup',filesep];
    end
    
    savename = [date,'_',num2str(latticeSize(1)),'x',num2str(latticeSize(2)),...
        '_U=',num2str(paramU),'_t=',num2str(paramT),'_nUp=',num2str(sum(basisStatesUp(1,:),2)),...
        '_nDown=',num2str(sum(basisStatesDown(1,:),2)),'PointGroup'];
    
    
    save([savedir,savename,'.mat'],'resultsPointGroup');
end

end



if flagDoNoSymmetry == 1% && nSites < 10
%% NO Symmetry case:
    tic
    
    HamNosym = HubbardHamiltonNoSymmetries(basisStatesUp, basisStatesDown, ...
        indNeighbors, paramU, paramT);
    
    E0_full = lanczos_myversion(HamNosym,nLanczosEigenvalues,lanczosPrecision,lanczosSteps);
    
    timeNoSymmetry = toc;
    
    resultsNoSymmetry = {'lattice size', 'U/t', '# spin UP', '# spin DOWN',...
        'time','E0','# states','PC info';  latticeSize, paramU/paramT, sum(basisStatesUp(1,:),2),...
        sum(basisStatesDown(1,:),2), timeNoSymmetry,E0_full,nStatesUp*nStatesDown,[cpuInfos,moreInfos]};%#ok

    
    if flagSaveResults == 1
        
        if isunix
            savedir = ['/afs/itp.tugraz.at/user/dobrautz/Dropbox/Master/Diplomarbeit/meineProgramme/Hubbard_ED/Results/noSymmetry',filesep];
        else
            savedir = ['D:\Program Files (x86)\Dropbox\Dropbox\Master\Diplomarbeit\meineProgramme\Hubbard_ED\Results\noSymmetry',filesep];
        end
        
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
    
end

clear all;
matlabpool close;
pctRunDeployedCleanup;

if isunix
    exit;
end

end

