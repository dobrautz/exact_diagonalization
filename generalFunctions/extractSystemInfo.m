function [systemFlag, nSites, paramU, nElec, representation, ...
    numWalkers, timeStep, twistedBCflag] = extractSystemInfo(loadDir)
% function to extract info about system like Hubbard U, lattice size etc.
% from the output string of writeFCIDUMP

% get name of file
dumpfileName = dir(fullfile(loadDir,'*FCIDUMP'));
% use only name of struct
dumpfileName = dumpfileName(1).name;

% if DUMPFILE is only called 'DUMPFILE' -> dont get any infos and abort.
if strcmp(dumpfileName,'FCIDUMP') == 1
    return
end

% if system is calculated for twiste boundary conditions an additional
% '_TBC_' specifier is in FCIDUMP name -> output it but remove this
if ~isempty(strfind(dumpfileName,'TBC_'))
    twistedBCflag = 1;
else
    twistedBCflag = 0;
end

dumpfileName = regexprep(dumpfileName, 'TBC_','');

% find underscore postions as those indicate differerent infos
barPositions = find(dumpfileName == '_');

% get system information
systemString = dumpfileName(barPositions(1):barPositions(2));

if ~isempty(strfind(systemString,'TBchain'))
    systemFlag = 'TB';
elseif ~isempty(strfind(systemString,'AIMchain'))
    systemFlag = 'AIM';
else
    systemFlag = 'Hubbard';
end

% get system size with formateed string read:
nSites = sscanf(systemString(2:end),'%i %*c %i');
% 'x' is always between x and y dim of lattice hence %*c to omit that char
nSites = nSites(1)*nSites(2);

% get parameter U value
paramUstring = dumpfileName(barPositions(2):barPositions(3));
paramU = sscanf(paramUstring(4:end),'%d');

% get number of electrons:
electronString = dumpfileName(barPositions(4):barPositions(5));
nElec = sscanf(electronString(8:end),'%i');

% get representation in which calculated
representationString = dumpfileName(barPositions(6):barPositions(7));
if ~isempty(strfind(representationString,'FALSE'))
    representation = 'RHF';
else
    representation = 'UHF';
end

% number of walkers
if ~isempty(strfind(dumpfileName,'numWalkers'))
    numWalkers = str2num(dumpfileName(barPositions(10)+1:barPositions(11)-1));
else
    numWalkers = nan;
end

% time step
if ~isempty(strfind(dumpfileName,'timeStep'))
    timeStep = str2num(dumpfileName(barPositions(12)+1:barPositions(13)-1));
else
    timeStep = nan;
end