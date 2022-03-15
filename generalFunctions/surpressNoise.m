function [ matrix ] = surpressNoise( matrix, surpressThreshold, flagRenorm)
%surpressNoise surpresses noise and renorm matrices (until now only for
%cubic matrices! and up to 4 dimensions!

if nargin < 3
    flagRenorm = 0;
end

matrixDim = ndims(matrix);

% surpress
% matrix(abs(matrix) < surpressThreshold) = 0;

matrix = real(matrix) .* (abs(real(matrix)) > surpressThreshold) + ...
    1i * imag(matrix) .* (abs(imag(matrix)) > surpressThreshold);

% and renorm 
if flagRenorm
    if matrixDim == 2
        
        matrix = bsxfun(@rdivide, matrix, sum(matrix.^2));
        
    elseif matrixDim == 3
        
        for iSize1 = 1:length(matrix)
            
            matrix(:,:,iSize1) = bsxfun(@rdivide, matrix(:,:,iSize1), ...
                sum(matrix(:,:,iSize1).^2));
            
        end
        
    elseif matrixDim == 4
        
        for iSize1 = 1:length(matrix)
            
            for iSize2 = 1:length(matrix)
                
                matrix(:,:,iSize1,iSize2) = bsxfun(@rdivide, matrix(:,:,iSize1,iSize2),...
                    sum(matrix(:,:,iSize1,iSize2).^2));
                
            end
            
            
        end
        
    end
end
