function [ hubbardPhases ] = HubbardPhase( replBasisStates , hopIndicationMatrix , hoppingIndicesX , hoppingIndicesY )
%HUBBARDPHASE calculates fermionic phases between 2 basis state matrices
%
% Input:     replBasisStates ...  augmented basis state set, copied for hopping
%           hopIndicationMatrix  matrix where 1s indicate hopping processes
%           hoppingIndicesX/Y..  site indizes where hopping occurs
%
% Output:    phases ...  phases from hubbard hopping
% 
% NOTE: 22.06.13: parallized loop over sites
%------------------------SVN Info------------------------------------------
% $Rev:: 82                                     $: Revision of last commit
% $Author:: dobrautz                            $: Author of last commit
% $Date:: 2013-06-22 19:51:27 +0200 (Sam, 22. J#$: Date of last commit
% -------------------------------------------------------------------------
%
% TODO: * this function is still bottleneck, try to improve it
%       * and same as in phase factor calculation of cubicSymmetries.m 
%           check if it is done right and make a better explanation how it 
%           works
%       * probably not vectorizeable ... only with extensive use of cell
%           arrays

% matrix which indicates sites where hopping occurs
hopIndicationMatrix = logical(hopIndicationMatrix)';

[nReplBasisStates,nSites] = size(replBasisStates); %numer of states and sites
hubbardPhases = zeros(nReplBasisStates,1); %container for phases

%preperation for general phase_fac function
% create indication matrix for flipped sites in transposed form
indSites = repmat((1:nSites)',1,nReplBasisStates); 

% possible way to do it without hoppingIndicesX/Y input: TODO check for
% speed:
% hoppingInd = reshape(indSites(hopIndicationMatrix),2,nReplBasisStates);
% hoppingInd = flipud(hoppingInd);
% hoppingInd = hoppingInd(:);

% combine indices in switched form, % have to sort it do correctly include 
% hopping over boundaries 
hoppingInd = [max(hoppingIndicesY,hoppingIndicesX),min(hoppingIndicesY,hoppingIndicesX)]'; 
hoppingInd = hoppingInd(:); %bring it in vector form

indSites(hopIndicationMatrix) = hoppingInd; % switch it

%% phase factor calculation
% bring it in normal form; need index transposed because of matlab first 
% row sorting:
indSites = flipud(indSites);

indReplBasisStates = (1:nReplBasisStates)';
% also need it in cell format for cellfun() application!
cellIndStates = num2cell(indReplBasisStates);

% NOTE: 01.06.13: did it!! borught it in vectorized form!!
for iSites = 1:nSites
    
    % matrix[iSites,nReplStates] containing vector 1:iSites as columns
    indUpToiSites = repmat((1:iSites)',1,nReplBasisStates);
    % pick elements up to(!!) element equal index: iSites from indSites and
    % cast it in right form:
    indUpToSiteIndex = reshape(indSites(indSites <= iSites),iSites,nReplBasisStates);
    % index up to which vector should go:
    indUpToVector = indUpToiSites(indUpToSiteIndex==iSites);
    % short form of this below:
    cellIndUpToVector = num2cell(indUpToVector);
    
    if iSites == 1
        cellIndUpToVector = cellIndUpToVector';
    end
    
    sumHoppedElectrons = cellfun(@(x,y) sum( replBasisStates(x,...
        indUpToSiteIndex(1:y-1,x)),2), cellIndStates , cellIndUpToVector);
    % see descripton below: this above only faster form!
%     % now i want elements from indUpToSiteIndex up to this value,
%     % but this depends for each state on its indSites value so this are
%     % different size vectors: solution use accumarray to cast it in a cell
%     % form
%     indexForUpToSiteIndex = accumarray(indReplBasisStates,indUpToVector,[],@(x) {1:(x-1)});
%     % equivalent to 1:(x(indUpToSiteIndex==iSites)-1) of old version below
%     % now want indUpToSiteIndex(indexForUpToSiteIndex) like below: use 
%     % cellfun()! but also need column index use indBasisStates in cell form
%     columnIndex = cellfun(@(x,y) {indUpToSiteIndex(x,y)}, indexForUpToSiteIndex, cellIndStates);
%     % now pick out values from basis with these indices use cellfun()
%     % because of needed cell format of columnIndex!
%     sumHoppedElectrons = cellfun(@(x,y) sum(replBasisStates(y,x),2), columnIndex, cellIndStates);
    % only counts if starting position had an 'electron/spin' there:
    checkIfElectron = replBasisStates(:,iSites) .* sumHoppedElectrons;
    % count exponential sum factors:
    hubbardPhases = hubbardPhases + checkIfElectron;
    % now have right indices now want to pick out corresponding values of
    % basis and sum over them
    
    %... TODO
    
end
% calc corresponding phase:
hubbardPhases = (-1).^hubbardPhases;

% this much easier version also works: just (-1) to the power of the sum of
% operators between the two hoppung positions! change this later
siteInd = 1:nSites;
phase = zeros(nReplBasisStates,1);
for i = 1:nReplBasisStates
    
    x = siteInd(hopIndicationMatrix(:,i));
    
    phase(i) = (-1)^sum(replBasisStates(i,(x(1)+1):(x(2)-1)));
    
end

assert(all(phase==hubbardPhases))
%-- NOTE: 01.06.13: vectorization done : this below is old version:--------
% vectorization probably not possible. and if only with extensive use of
% cell arrays. 
% secon idea to use a unique comman beforehand and only do it for unique
% elements test: 
% NOTE: 19.05.13 not useful too, because not only unique indices are
% important but combination of unique indices and binary basis
% representation! and thats seldom! 

% indSites = indSites';
% indSites = fliplr(indSites);

% for iReplStates = 1:nReplBasisStates  % loop over all states
%     
%     % pick out specific indices and states
%     specIndSites = indSites(iReplStates,:);
%     specBasisStates = replBasisStates(iReplStates,:);
%     
%     % and loop over sites, similar to my function phase_fac look there for 
%     % explanation, or now in cubicSymmetries.m since rework 
%     for iSites = 1:nSites 
%         
%         indUpToSiteIndex = specIndSites(specIndSites <= iSites);
%         x = 1:iSites;
%         y = indUpToSiteIndex(1:(x(indUpToSiteIndex==iSites)-1));
%         z = sum(specBasisStates(:,y),2);
%         n = specBasisStates(:,iSites) .* z;
%         
%         hubbardPhases(iReplStates) = hubbardPhases(iReplStates) + n;
%         
%     end
%     
% end
% 
% hubbardPhases = (-1).^hubbardPhases;


