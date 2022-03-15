function [ cumulIndex ] = calcCumulativeIndex( normHubbardStates, nBasisStatesDown )
%CALCCUMULATIVEINDEX calculates from referencing style save method of
% normHubbardStates the cumulative Index

% determine which states just reference to whole basis have a -1 as only
% element:
% creates cell array ( from option: 'uniformoutput'=false) indicating if
% elements are negative. 
indNegative = cellfun(@(z) z <0,normHubbardStates,'uniformoutput',0);
% with all function: gives one TRUE if all elements are TRUE:
indNegative = cellfun(@all ,indNegative); 

% determine length of all elements
elementLength = cellfun(@length,normHubbardStates);
% first elements as long as whole basis
firstElementWholeBasis = elementLength(elementLength == nBasisStatesDown);
if ~isempty(firstElementWholeBasis) 
firstElementWholeBasis = firstElementWholeBasis(1);
end

% replace values of referencing states:
elementLength(indNegative) = firstElementWholeBasis;
% calc cumulative sum:
cumulIndex = [0;cumsum(elementLength)];



end

