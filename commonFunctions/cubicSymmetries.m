function [ basisStates , symmetryPhases ] = cubicSymmetries( basisStates , ...
    symmetryValue , latticeSize , varargin)
%cubicSymmetries Combination of former mirror... and translation...
% functions, IMPORTANT the inputname of the second input decides which
% symmetry is applied. e.g. py is a mirror symmetry on Y axis , rx is a
% X-axis translation etc. , except on explicetly calls the function with
% a third input argument, which is a string definig the to applied symmetry
%
% at least 2 input arguments are needed, if no 3rd string input argument 
% defining the symmetry operation is provided, the inputname(2) of the 
% second inpput is taken. this means the function has to be called in the 
% right way!!! IMPORTANT!!!
% 
% Input:  basisStates       binary form basis on which symmetry is applied
%         symmetryValue     symmetry Eigenvalue !! inoutname of this decides
%                           which symmetry operation is applied:
%                           px      mirror X axis
%                           py      mirror Y axis
%                           rx     translation X axis etc. ..
%                           a zero value here causes nothiun to happen!
%         latticeSize       for rectangular lattices i now need also
%                           lattice as input 
%         string 4th input  of form 'px','py' etc. see above. for cases 
%                           where function cant be called in right form
% 
% Output: newBasisStates    applied symmetry on basis
%         symmetryPhases    correspondind fermionic phases
% 
% NOTE: 25.06.13: change to also include rectangular lattices-> need
% another input of lattice size
% finished changes -> unit_testing todo, and take care in translation cases
% other as matlab usual -> latticeSize = [Lx;Ly] -> reshape(x,Ly,Lx) !! 
% NOTE: 25.06.13: switched px and py now makes more sense in combination
% with pd and pe cases!
% NOTE: 10.07.13: change latticeSize = [nSitesY;nSitesX] like matlab row
% column style
% 
%------------------------SVN Info------------------------------------------
% $Rev:: 60                             $: Revision der letzten Übertragung
% $Author:: dobrautz                    $: Autor der letzten Übertragung
% $Date:: 2013-06-10 11:39:59 +0200 (Mo#$: Datum der letzten Übertragung
% -------------------------------------------------------------------------
% 
% TODO: * check if JIT-Accelerator is running (29.04.2013)
%       * make documentation of phase factor calc

%MIRROR.. applies the respective symmetry operation to binary
%heisenberg/hubbard basis and calculates the respective fermionic sign
%use py = 0 to surpress symmetry application!

[nBasisStates,nSites] = size(basisStates); %size of basis
symmetryPhases = zeros(nBasisStates,1);   % cont. var. for phase fac
% 1D lengths NOTE: 25.06.13: change to rectangular lattices 
nSitesY = latticeSize(1);

%% identify which symmetry operation is to be applied
if numel(varargin) > 0
    
    input_name = varargin{1};
    
else
    
    input_name = inputname(2);
    
end

if symmetryValue == 0% dont do anything if py = 0 -> sign = 1
    
    symmetryPhases = 1.^symmetryPhases;
    
else    %else do something
    %% Symmetry Operation
    % look up in documentation for derivations of respective forms below:
    % apply chosen symmetry operation
    switch input_name
        
        case 'py' %Mirror Y-Axis
            %--------------------------------------------------------------------
            
            z = nSitesY:-1:1;
            z1 = 0:nSitesY:nSites-nSitesY;
            % bsx version:
            switchIndex = bsxfun(@plus,z,z1')'; 
            switchIndex = switchIndex(:)';
%             z = repmat(z,1,nSitesX);
%             z1 = repmat(z1,nSitesY,1);
%             z1 = z1(:)';
%             switchIndex = z + z1;

            %--------------------------------------------------------------------
            
        case 'px' %Mirror X-Axis
            %---------------------------------------------------------------
            
            z = 1:nSitesY;
            z1 = nSites-nSitesY:-nSitesY:0;
            % bsx version:
            switchIndex = bsxfun(@plus,z,z1')'; 
            switchIndex = switchIndex(:)';
%             z = repmat(z,1,nSites1D);
%             z1 = repmat(z1,nSites1D,1);
%             z1 = z1(:)';
%             switchIndex = z + z1;
            %---------------------------------------------------------------
        
        case 'pe' %Mirror E-Axis
            
            % -------------------------------------------------------------
            z = nSites:-nSitesY:nSitesY;
            z1 = 0:nSitesY-1;
            % bsx version:
            switchIndex = bsxfun(@minus,z,z1')'; 
            switchIndex = switchIndex(:)';
%             z = repmat(z,1,nSites1D);
%             z1 = repmat(z1,nSites1D,1);
%             z1 = z1(:)';
%             switchIndex = z - z1;
            % -------------------------------------------------------------    
        case 'pd' %Mirror D-Axis
            
            %-------------------------------------------------------------
            switchIndex = reshape(1:nSites,nSitesY,nSitesY);
            switchIndex = switchIndex.';
            switchIndex = switchIndex(:).';
            % -----------------------------------------------------------
          
        case 'rx' %Translation X-Axis
            nSitesX = latticeSize(2);
            %---------------------------------------------------------------------
            switchIndex = reshape(1:nSites,nSitesY,nSitesX);
            switchIndex = circshift(switchIndex,[0,symmetryValue]);
            switchIndex = switchIndex(:)';
            %---------------------------------------------------------------------
            
        case 'ry' %Translation Y-Axis
            
            nSitesX = latticeSize(2);
            %------------------------------------------------------------------
            switchIndex = reshape(1:nSites,nSitesY,nSitesX);
            switchIndex = circshift(switchIndex,symmetryValue);
            switchIndex = switchIndex(:)';
            %------------------------------------------------------------------
            
        case 'rx1D' % Translation X axis 1 dimensional
            
            switchIndex = 1:nSites;
            switchIndex = circshift(switchIndex,[0,symmetryValue]);
            
    end
    
    % copy old basis
    basisStatesCopy = basisStates;
    % apply symmetry to basis
    basisStates(:,1:nSites) = basisStates(:,switchIndex);
    
    %% Phase factors
    
    % found one function which works for all symmetry operations
    % maybe i can include this in the hamiltonian construction to, since its
    % similar
    % derivation for this in documentation
    flippedSwitchIndex = fliplr(switchIndex);
    
    parfor iSites = 1:nSites %loop over sites
        
        indSwitchIndSmallerSiteInd = flippedSwitchIndex(flippedSwitchIndex <= iSites); %#ok
        indUpTo = 1:numel(indSwitchIndSmallerSiteInd);
        
        y = indSwitchIndSmallerSiteInd(1:indUpTo(indSwitchIndSmallerSiteInd==iSites)-1);
        x = sum(basisStatesCopy(:,y),2); %#ok
        n = basisStatesCopy(:,iSites) .* x;
        symmetryPhases = symmetryPhases + n;
        
    end
    
    % phase 
    symmetryPhases = (-1).^symmetryPhases;

end
