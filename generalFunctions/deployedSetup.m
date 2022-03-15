function [saveDir, inputData] = deployedSetup
%deployedSetup loads input data and creates save dir in deployed tag
%releases of function and closes existing matlab pool sessions and sets
%nohup execution flags
% 
% Output:   inputData       
%           savedir
% 
% Notes: 
% NOTE: 10.11.13: created
%                 creates results folder in folder of called function in
%                 this way easier to setup generally for all deployed
%                 functions

% has to called directly by function for dbstack(1) to work
callingFile = dbstack('-completenames',1);
currentFolder = regexprep(callingFile.file,[callingFile.name,'.m'],'');

% currentFolder = regexprep(mfilename('fullpath'),mfilename,'');
inputFileName = dir([currentFolder,'*input*.txt']);
saveDir = [currentFolder,'Results',filesep];

% only load inputData if requested or else probably no input file
if nargout > 1
    inputData = load([currentFolder,inputFileName.name]);
end

if ~exist(saveDir,'dir')
    mkdir(saveDir);
end

% on unix systems for nohup execution call this function right here
if isunix && (isdeployed || ismcc)
    system('export MCR_CACHE_ROOT=/tmp/dobrautz')
end

% check if matlab pool is open and close if any are open
if matlabpool('size') > 0
    matlabpool close;
end
