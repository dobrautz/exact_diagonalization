function [whichSymmetry, symmetrySector] = getSymmetryString( symmetrySector,lattice )
%GETSYMMETRYSTRING determines symmetry defining string dependent on
%symmetry.pointGroup field array. CAUTION changes symmetrySector struct!!

if all(symmetrySector.pointGroup ~=0)
    
    assert(lattice.size(1)==lattice.size(2))
    
    whichSymmetry = '';
    
    
elseif numel(symmetrySector.pointGroup ~= 0) == 2
    
    whichSymmetry = '';
    
    symmetrySector.pointGroup = symmetrySector.pointGroup(symmetrySector.pointGroup~=0);
    
elseif symmetrySector.pointGroup(1) ~= 0;
    
    assert(lattice.size(1)==lattice.size(2))
    
    whichSymmetry = 'pd';
    
    symmetrySector.pointGroup = symmetrySector.pointGroup(symmetrySector.pointGroup~=0);
    
elseif symmetrySector.pointGroup(2) ~= 0;
    
    whichSymmetry = 'py';
    
    symmetrySector.pointGroup = symmetrySector.pointGroup(symmetrySector.pointGroup~=0);
    
elseif symmetrySector.pointGroup(3) ~= 0;
    
    whichSymmetry = 'px';
    
    symmetrySector.pointGroup = symmetrySector.pointGroup(symmetrySector.pointGroup~=0);
    
end


