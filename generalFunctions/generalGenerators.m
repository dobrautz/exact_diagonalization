function E = generalGenerators(basicGenerator,xIndex,yIndex)
% GENERALGENERATORS recursively calculates a general unitary group
% generator E_ij through usage of the commutator relations

if (yIndex - xIndex) == 1
    
    E = basicGenerator{xIndex};
    
elseif (yIndex - xIndex) == 2 
    
    E = basicGenerator{xIndex}*basicGenerator{xIndex+1} - ...
        basicGenerator{xIndex+1}*basicGenerator{xIndex};
    
else 
    
    X = generalGenerators(basicGenerator, xIndex, yIndex-1);
    
    E = X*basicGenerator{yIndex-1} - basicGenerator{yIndex-1}*X;
    
end
