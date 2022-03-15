function H_down = hamiltonDownNoSyms( basisStatesDown , nBasisStatesUp , paramT , indNeighbors )
%HAMILTONDOWNNOSYMS calculates down spin hopping part of hamiltonian
%without symmetries
% NOTE: 30.05.13:  make it more beatiful

% number of states, sites and dimension
[nBasisStates,nSites] = size(basisStatesDown);
[nDims,~] = size(indNeighbors);
% for binary2dez conversion
bin2dez = nSites-1:-1:0;   
bin2dez = (2.^bin2dez)';  

%container for matrix indices (need cell structure because of 2 dimensions)
[xIndHubbard,yIndHubbard,phasesHubbard] = deal(cell(1,nDims)); 

%-------------------------- note 10.04.2013-------------------------------
% very confusing variable names and structure of this program, write nice
% version of this sometimes
%--------------------------old notes---------------------------------------
% this part below calculates intern index distribution for each UP Spin
% repr., since they are ordered: UP: MSB, DOWN:LSB. and for non vanishing
% overlap of two DOWN spin parts the UP Spin parts have to be the same
% and since (for now atleast, maybe due to change 11.04.2013) there is
% always a comp. repr. created if one just switches the DOWN Spin part of a
% state, the new states are easy to handle.
%--------------------------------------------------------------------------
% program part below is written very storage effective for huge vectors and
% matrices, see version: heisenberg_hamiltonian_pretty for nice
% documentation, TODO

% parfor (iDims = 1:nDims, 1) %loop over Dimensions
for iDims = 1:nDims % dont forget clear command below
    
B2 = basisStatesDown(:,indNeighbors(iDims,:));      %binary form of neighbors in DIM i
B2 = xor(basisStatesDown,B2);         % determines sites where electrons can hop

%gets indices of sites where electrons can hop, but:
[f,d] = find(B2);         %f zeilen, d spalten ,nach spalten sortiert oO
[xIndHubbard{iDims},f] = sort(f);   % x_ind waere zeilen index von H_heisenberg

d = d(f);                 %sort d same way as f, is first site index of 2 sites where hopping occurs
f = indNeighbors(iDims,d)';             %index of neighbour of to be flipped site
f1 = [d;f];               %combined index for flipped sites
d1 = [1:sum(B2(:)),1:sum(B2(:))]';%two times index from 1:N for indexing the flipped sites
                                  %since i create in next row a new matrix
                                  %of binary form of, where basis is
                                  %reproduced so many times as hops can
                                  %occur

F = basisStatesDown(sum(B2,2)==0,:)*bin2dez;           %fehlende zahlen-indizes: as some hops to integer values are not allowed
                                    %but one still needs this values for
                                    %correct ordering of new states to
                                    %determine y_ind value
                                    
B2 = basisStatesDown(xIndHubbard{iDims},:);     %basis so oft reproduziert wie gswitcht wird, see above

d1 = sparse(d1,f1,1);              %indication matrix for switches

% now calculate the fermionic phases of hamilton hops:
[phasesHubbard{iDims}] = HubbardPhase(B2,d1,d,f);%maybe have to sort this if x and y are sorted
%TO DO hubbardphase is bottleneck!! make it better!

d = xor(d1,B2);                  %new states, thats enough

d = d*bin2dez;                       %corresponding integer values
f = [d;F];                      %include missing integer values for comparison

[~,~,B2] = unique(f);           %check them! thats enough! since standard call of unique
                                % sorts the unique-values, indices B2 give
                                % me right positions of new created states 
yIndHubbard{iDims} = B2(1:numel(xIndHubbard{iDims}));       %only include the one you need, possible since missing integer
                                        % values were included at the end
                                        % of list
                         
clear f B2 d d1 f1 

end


%Indizes der nichtdiagonalelemente
xIndHubbard = cell2mat(xIndHubbard(:));  
yIndHubbard = cell2mat(yIndHubbard(:)); 
phasesHubbard = cell2mat(phasesHubbard(:)); 

%--------------------------old notes--------------------------------------
% fuer jeden rep. gibts nun diese interne indizes verteilung der down spins
% this means one has to replicate the indices for each repr. while
% increasing their values by the number of DOWN spin basis states for each
% repr. 
% x is the starting state index, y the end state index
% still storage effective code, since vectors become huge!
%--------------------------------------------------------------------------
% ------------------------------ old way:----------------------------------
% % increasing vector for each repr. replicated in right form for each new state
% indRepr = 0:nBasisStatesUp-1;
% indRepr = repmat(indRepr,numel(xIndHubbard),1);   
% indRepr = indRepr(:);
% % replicate x_ind for each repr. +calc x in x = i*N + k; confusing var.names! 
% xIndHubbard2 = repmat(xIndHubbard,nBasisStatesUp,1);     
% % vector with # of down states times onld indices
% indRepr = nBasisStates*ones(numel(xIndHubbard2),1) .* indRepr;    
% 
% % calc x in x = i*N + k; confusing var. names!
% xIndHubbard2 = xIndHubbard2 + indRepr; 
% 
% % same with y_ind
% yIndHubbard2 = repmat(yIndHubbard, nBasisStatesUp,1);    
% yIndHubbard2 = yIndHubbard2 + indRepr;
%------------- NOTE: 31.05.13: new way of doing it with bsxfun()!---------
% increasing vector for each repr. replicated in right form for each new state
indRepr = (0:nBasisStatesUp-1) * nBasisStates;
% thats enough: bsxfun does the right repamt functions of above!
xIndHubbard = bsxfun(@plus, xIndHubbard, indRepr);
xIndHubbard = xIndHubbard(:);
yIndHubbard = bsxfun(@plus, yIndHubbard, indRepr);
yIndHubbard = yIndHubbard(:);

%------------------------note 11.04.2013-----------------------------------
% as of this date the hubbard norm is a huge vector containing the norm
% values of each hubbard cycle repr. in order form for UP spin repr.
%--------------------------------------------------------------------------

% number of hubbard basis states
nHubbardStates = nBasisStatesUp*nBasisStates;

% replicate phases (are same for every repr.) and multiply with t to get H
H_down = -paramT * repmat(phasesHubbard,nBasisStatesUp,1); 

% pseudo values for sparse creation, assures that matrices are of right size
[xIndHubbard(end+1),yIndHubbard(end+1)] = deal(nHubbardStates);

clearvars -except *IndHubbard H_down

H_down = sparse(xIndHubbard,yIndHubbard,[H_down;0]);

clear *IndHubbard

% assure hermiticity
H_down = (H_down + H_down')/2;