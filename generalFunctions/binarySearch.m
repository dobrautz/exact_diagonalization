function [index, flagNotFound] = binarySearch(list,key)
%binSearch is a simple binary search function on a sorted list of integers
% 
% Input : list          vector of sorted unique integes
%         key           value to be found
% 
% Output: index         index of element if found, if not found gives index
%                       of last compared element and the flag indicates if
%                       its one lower or higher
%         flagNotFound  flag indicating if element was in list:
%                       0.. element was in list
%                       1.. element is bigger than last compared one
%                       -1.. elements is smaller than last compared one
% 
% NOTE: 14.10.13: created: source:jagger.berkeley.edu/~pack/e77/Searching.ppt?

left = 1;
right = length(list);

while right > left
    
    mid = floor((left+right)/2);
    
    if list(mid) < key
        
        left = mid+1;
        
    else
        
        right = mid;
    end
end

index = left;

if list(left)==key
    
    flagNotFound = 0;
    
elseif list(left) < key
    
    flagNotFound = 1;
    
elseif list(left) > key
    
    flagNotFound = -1;
    
end
