function results = tableUpdate(results,E0,V0,kx,ky,varargin)
%TABLEUPDATE Updates the results-table with results from lanczos
% diagonalisation for specific translational and point symmetries
%
%------------------------SVN Info------------------------------------------
% $Rev:: 64                             $: Revision der letzten Übertragung
% $Author:: dobrautz                    $: Autor der letzten Übertragung
% $Date:: 2013-06-13 17:46:13 +0200 (Do#$: Datum der letzten Übertragung
% -------------------------------------------------------------------------



%container
symmetryValues = zeros(1,6);

symmetryValues(1) = kx; symmetryValues(2) = ky;

% check for input:
if numel(varargin) == 3
    
    symmetryValues(3) = varargin{1};
    symmetryValues(4) = varargin{2};
    symmetryValues(5) = varargin{3};
    
elseif numel(varargin) == 2
    
    symmetryValues(3) = varargin{1};
    symmetryValues(4) = varargin{2};
    
elseif numel(varargin) == 1
    
    input_name = inputname(6);
    input_value = varargin{1};
    if strcmp(input_name,'px');
        
        symmetryValues(3) = input_value;
        
    elseif strcmp(input_name,'py');
        
        symmetryValues(4) = input_value;
        
    elseif strcmp(input_name,'pd');
        
        symmetryValues(5) = input_value;
        
    else
        
        symmetryValues(6) = input_value;
    end
    
end

if isempty(results)==1 %bei erstem Aufruf
    
    results = cell(1,8);
    
    
    for i = 1:6
        results{i} = symmetryValues(i);
    end
    results{7} = E0;
    results{8} = V0;
    
else
    
    [j,~] = size(results);
    j = j+1;

    for i = 1:6
        results{j,i} = symmetryValues(i);
    end
    results{j,7} = E0;
    results{j,8} = V0;
    
    
    
    
end

