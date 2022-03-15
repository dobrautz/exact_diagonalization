function [ FCIdata, FCIcomments ] = readFCIMCStats( fileLocation, varargin )
%READFCIMCSTATS reads in FCIMCStats files without possible columns
% containing 'Infinity' due to Matlab being incompatible with it.
% 2 possible implementations using fscanf and textscan
% fscanf version only for real NECI output/not yet kNECI compatible
% textscan output is cell!!
FCIMCStats = fopen(fileLocation,'r');

% first 2 lines are comments
FCIcomments{1} = fgetl(FCIMCStats);
FCIcomments{2} = fgetl(FCIMCStats);

% third line is dummy line
fgetl(FCIMCStats);

if nargin < 2 || strcmp(varargin{1},'fscanf')
    % fscanf implementation (standard):
    % define input read-string with skipped lines 
    % 11.Proj.E.ThisCyc
    % 21.HFInstShift
    % 23.Tot-Proj.E.ThisCyc
    % 30.AbsProjE
    FCIMCStatsString = '%i %f %i %f %i %i %i %i %f %f %*s %i %i %e %i %e %f %i %e %f %*s %f %*s %f %f %f %f %f %f %*s %i';
    FCIdata = fscanf(FCIMCStats, FCIMCStatsString, [27,inf]);
    FCIdata = FCIdata';

else 
    % textscan implementation
    FCIMCStatsString = [repmat('%n ',1,10),'%s ',repmat('%n ',1,11),'%s ',repmat('%n ',1,8)];

    FCIdata = textscan(FCIMCStats,FCIMCStatsString,'treatAsEmpty',{'Infinity','-Infinity','Infi','-Infi'},'delimiter',' ','commentStyle', '#','MultipleDelimsAsOne',1);

end

fclose(FCIMCStats);
