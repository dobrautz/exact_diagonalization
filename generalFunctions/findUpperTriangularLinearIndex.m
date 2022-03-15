function linearIndex = findUpperTriangularLinearIndex(input)
% finds the matrix size linearIndex*linearIndex if only the number of
% elements in the upper triangular parts are given.

linearIndex = 1;
nElements = 0;

while linearIndex < 100 % bigger matrices not expected to be needed
    nElements = nElements + linearIndex;
    
    if input == nElements
        break
    else
        linearIndex = linearIndex + 1;
    end
end

        