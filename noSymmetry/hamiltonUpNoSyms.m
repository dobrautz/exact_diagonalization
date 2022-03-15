function H_up = hamiltonUpNoSyms( basisStatesUp, basisStatesDown , paramT , indNeighbors )
%HAMILTONUPNOSYMS calculates the up spin hopping part of hamiltonian witout
%applied symmetry operations

%number of spin down states, sites and up spin repr.
[nBasisStatesDown,nSites] = size(basisStatesDown);  
[nBasisStatesUp,~] = size(basisStatesUp);    

[nDims,~] = size(indNeighbors);

bin2dez = nSites-1:-1:0;   %for binary2dez conversion
bin2dez = (2.^bin2dez)';  

%container for matrix indices need cell structure because of 2 dimensions
[xIndHubbard,yIndHubbard,phasesHubbard] = deal(cell(1,nDims));

%-------------------------- note 10.04.2013-------------------------------
% very confusing variable names and structure of this program, write nice
% version of this sometimes
%-------------------------------------------------------------------------

% program part below is written very storage effective for huge vectors and
% matrices, see version: heisenberg_hamiltonian_pretty for nice
% documentation, TO DO

%-------------------------note 11.04.2013----------------------------------
% as opposed to the spin DOWN part of the hamilton matrix, one needs to
% take more care in calculation of the UP spin part. first of all: there is
% such an easy index distribution, since the UP spin part is the MSB part
% of the binary form.( maybe with intelligen switches between MSB and LSB
% there is some kind of order which can be used, not done yet) 
% and secondly: now there exist the possibility that a non comp. repr. or
% no repr. at all is created by a hamilton electron hop. so one has t otake
% greater care in calculations of state overlaps. see notes below
%--------------------------------------------------------------------------
% keyboard

% parfor(iDims = 1:nDims, 1) %loop over Dimensions
for iDims = 1:nDims  %dont forget clear command!!

B2 = basisStatesUp(:,indNeighbors(iDims,:)); % binary form of neighbors in DIM i
B2 = xor(basisStatesUp,B2);              % determines sites where electrons can hop
% gets indices of sites where electrons can hop, but:
[f,d] = find(B2);   %f zeilen, g spalten ,nach spalten sortiert oO

% for UP Spin part one has to check if there is a repr. at all! or a repr.
% with no hopping possibilities in a specific dimension
if ~isempty(f)      %if there is hopping possible...
    
% and a stupid matlab inconsitency has to be checked:  
% wenn nur eine zeile in B_repr macht mir find() einen zeilenvektor!!! wtf
if nBasisStatesUp == 1    
    f = f'; d = d'; % make column vector of it
end

[xIndHubbard{iDims},f] = sort(f);     % x_ind waere zeilen index von H_heisenberg

% sort d same way as f, is first site index of 2 sites where hopping occurs
d = d(f);            
f = indNeighbors(iDims,d)';          %i ndex of neighbour of flipped site

f1 = [d;f];                       % combined index for flipped sites
d1 = [1:sum(B2(:)),1:sum(B2(:))]';% two times index from 1:N for indexing the flipped sites
                                  % since i create in next row a new matrix
                                  % of binary form of, where basis is
                                  % reproduced so many times as hops can
                                  % occur
                                  
F = basisStatesUp(sum(B2,2)==0,:)*bin2dez;  %fehlende zahlen-indizes: as some 
                                        %hops to integer values are not allowed
                                        %but one still needs this values for
                                        %correct ordering of new states to
                                        %determine y_ind value
                                    
B2 = basisStatesUp(xIndHubbard{iDims},:);  %basis so oft reproduziert wie gswitcht wird 

d1 = sparse(d1,f1,1);         %indication matrix for switches

% now calculate the fermionic phases of hamilton hops:
% calc. hopping phases, maybe have to sort this if x and y are sorted
[phasesHubbard{iDims}] = HubbardPhase(B2,d1,d,f); 
% TODO hubbardphase is bottleneck!! make it better!

d = xor(d1,B2);                  % neue Zustaende!

d = d*bin2dez;   % corresponding integer values
f = [d;F];   % include missing integer values for comparison

[~,~,B2] = unique(f);% check them! thats enough!since standard call of unique
                     % sorts the unique-values, indices B2 give
                     % me right positions of new created states 
                                
yIndHubbard{iDims} = B2(1:numel(xIndHubbard{iDims}));%only include the one you need, possible since missing integer
                                 % values were included at the end of list

%---------------------note 11.04.2013--------------------------------------
% the x indices are now in the range of the number of repr. but the y
% indices range til the number of down spin basis states!
%--------------------------------------------------------------------------

clear f B2 d d1 f1
 
end

end

% Indizes des nichtdiagonalelemente
xIndHubbard =cell2mat(xIndHubbard(:));  
yIndHubbard = cell2mat(yIndHubbard(:)); 
phasesHubbard = cell2mat(phasesHubbard(:)); 

% and now have to calc. index x = i*N + k and phase and norm in right
% combination: + have to take care of 2D source of x and y indizes!!! TO DO

%--------------------note 11.04.2013---------------------------------------
% similar as in 'HubbardHamDown' the now found x,y indices has to be
% transformed for the combined hubbard basis states. but there are major
% differences between the two. 
% first: the UP spin repr. part of the basis is the MSB part of the binary
% representation. so they have 

%--------------------------note 12.04.2013---------------------------------
% das unten besser erklaeren:

% # of hamilton indices
% nNewStates = numel(xIndHubbard);       

% k index; also need this vector to get k' values since sort does it in a 
% stupid way below in loop, faster way possible TODO
indBasis = 1:nBasisStatesDown;        %index for down basis states
% ------------------------------ old way:----------------------------------
% %repl this as often as there are indices for a down spin state
% matIndBasis = repmat(indBasis,nNewStates,1); 
% 
% % repl x index for each down spin state;
% xIndHubbard = repmat((xIndHubbard-1)*nBasisStatesDown+1,1,nBasisStatesDown); 
% % operation (x_ind-1)*nr_basis+1 because... TODO erklaeren
% 
% 
% xIndHubbard = (xIndHubbard + matIndBasis)'-1; 
% xIndHubbard = xIndHubbard(:); %get x = i*N+K
% 
% yIndHubbard = repmat((yIndHubbard-1)*nBasisStatesDown+1,1,nBasisStatesDown);
% yIndHubbard = (yIndHubbard + matIndBasis)'-1; 
% yIndHubbard = yIndHubbard(:); %get x = i*N+K
%------------- NOTE: 31.05.13: new way of doing it with bsxfun()!---------
xIndHubbard = (xIndHubbard-1) * nBasisStatesDown + 1;
xIndHubbard = bsxfun(@plus, xIndHubbard, indBasis)' - 1;
xIndHubbard = xIndHubbard(:);

yIndHubbard = (yIndHubbard-1) * nBasisStatesDown + 1;
yIndHubbard = bsxfun(@plus, yIndHubbard, indBasis)' - 1;
yIndHubbard = yIndHubbard(:);
% TODO: erklaer das besser
%--------------------------------------------------------------------------

% clear matInd

%now calc old norm and phase in right combination:

%---------note 12.04.2013--------------------------------------------------
% is this the right form for the phases??? yes it is!:
% phases are sorted in right way for each UP spin repr. i have one of this 
% phases nor for each DOWN spin part of basis, 
% --------------------------------------------------------------------------

% replicate phases (are same for every repr.)
phasesHubbard = repmat(phasesHubbard',nBasisStatesDown,1); 
phasesHubbard = phasesHubbard(:);
% number of hubbard basis states
nHubbardStates = nBasisStatesUp*nBasisStatesDown;

H_up = -paramT .* phasesHubbard; 

% pseudo values for sparse creation, assures that matrices are of right size
[xIndHubbard(end+1),yIndHubbard(end+1)] = deal(nHubbardStates); 

clearvars -except *IndHubbard H_up

H_up = sparse(xIndHubbard,yIndHubbard,[H_up;0]);

clear *IndHubbard

% assure hermiticity
H_up = (H_up + H_up')/2;
