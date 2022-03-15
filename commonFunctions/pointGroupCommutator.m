function [ flagTrue ] = pointGroupCommutator( varargin )
%POINTGROUPCOMMUTATOR tests point group commutation relations
% NOTE: 25.06.13: doesnt work anymore with change of cubicSymmetries to
% include rectangular lattices
x = 1:9;

if nargin == 2 
    
    string1stSymmetry = varargin{1};
    string2ndSymmetry = varargin{2};
    
    flagTrue = all(cubicSymmetries(cubicSymmetries(x,1,string1stSymmetry),1,string2ndSymmetry) == ...
        cubicSymmetries(cubicSymmetries(x,1,string2ndSymmetry),1,string1stSymmetry));
    
else
    if nargin == 4;
        transValue = 1;
    else
    
        transValue = varargin{5};
    end
    
    string1stSymLeft = varargin{1};
    string2ndSymLeft = varargin{2};
    
    string1stSymRight = varargin{3};
    string2ndSymRight = varargin{4};
    
    flagTrue = all(cubicSymmetries(cubicSymmetries(x,1,string1stSymLeft),1,string2ndSymLeft)== ...
        cubicSymmetries(cubicSymmetries(x,1,string1stSymRight),transValue,string2ndSymRight));
    
end

