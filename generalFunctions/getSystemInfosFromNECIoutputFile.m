function [GSenergy, nOrbitals, nElectrons, refEnergy, currentProjEnergy] = ...
    getSystemInfosFromNECIoutputFile(directory, isUzeroFlag)
% reads out system information from a NECI output file in 'directory'


% determine if standard called output files are here
for i = 1:3
    switch i
        case 1
            % since 3 standard output file names(normal, condor, dcluster)
            outputfile = dir(fullfile(directory,'OUT.*'));
            outputName  = pickOldestOutputfile( outputfile );
            if ischar(outputName)
                break;
            end
        case 2
            outputfile = dir(fullfile(directory,'*.out'));
            outputName  = pickOldestOutputfile( outputfile );
            if ischar(outputName)
                break;
            end
        case 3
            outputfile = dir(fullfile(directory,'*.o*'));
            outputName  = pickOldestOutputfile( outputfile );
            if ischar(outputName)
                break;
            end
    end
end

assert(ischar(outputName),'no output file automatically found...')

if isUzeroFlag == 1
    % for U=0 case sometimes NECI crashes when writing a POPSFILE so take
    % other quantitiy in this case
    disp('assuming U = 0 taking <D0|H|D0> as GS energy!')
    [~,GSenergy] = system(['grep " <D0|H|D0>= " ',fullfile(directory,outputName)]);
    GSenergy = sscanf(GSenergy,' <D0|H|D0>=  %f');
    
    GSenergy = GSenergy(1);
    refEnergy = GSenergy;
else
    
    [~,GSenergy] = system(['grep " Summed approx E(Beta)= " ',fullfile(directory,outputName)]);
    
    GSenergy = sscanf(GSenergy,' Summed approx E(Beta)=  %f');
    
    [~,refEnergy] = system(['grep "Current reference energy" ',fullfile(directory,outputName)]);
    refEnergy = sscanf(refEnergy,' Current reference energy %f');

end

[~, nOrbitals] = system(['grep "  NUMBER OF SPIN ORBITALS IN BASIS : " ', fullfile(directory,outputName)]);
nOrbitals = sscanf(nOrbitals(38:end),'%i')/2; % half since # of spin-orbitals are in OUT.1

[~, nElec] = system(['grep "  NUMBER OF ELECTRONS :" ', fullfile(directory,outputName)]);
nElectrons = sscanf(nElec(25:end),'%i');

% get also current projected Energy for each starting vector
% write grep output in tempFile
system(['grep -B 3  "Writing a 64-bit POPSFILE..." ', fullfile(directory,outputName),' > ',...
    fullfile(directory,'tempFile')]);
% read back in tempFile line per line until EOF
fileID = fopen(fullfile(directory,'tempFile'),'r');
count = 0;

while true
    string = fgetl(fileID);
    % filter out non data entries
    if isempty(string)
        continue
    elseif strcmp(string,'*********************************')
        continue
    elseif strcmp(string,'Writing a 64-bit POPSFILE...')
        continue
    elseif strcmp(string,'--')
        continue
    elseif ~ischar(string)
        break
    elseif strcmp(string(1:11),' Totwalkers')
        continue
    else
        data = sscanf(string,'%i %e %i %f %i %i %i %i %e %e %e %i %i %e %i %e');
        count = count + 1;
        currentProjEnergy(count) = data(11);%#ok
    end
end
fclose(fileID);

delete(fullfile(directory,'tempFile'));

