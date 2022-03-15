function renameMatlabWindow
%RENAMEMATLABWINDOW Renames Matlab window title
% only do it at uni linux for now:
if isunix
   
    % hostname:
    [~,host] = system('hostname');
    % number of cores
    numCores = num2str(feature('numCores'));
    % cpu clock
    [~,cpuClock] = system('cat /proc/cpuinfo');
    cpuClock = cpuClock(find(cpuClock=='@',1):find(cpuClock=='@',1)+8);
    % size of physical memory:
    [~,memSize] = system('cat /proc/meminfo');
    memSize = num2str(round(str2num(regexprep(memSize(find(memSize==':',1)+1:...
        find(memSize=='k',1)-1),' ',''))/1000000));
    % set display name
    jDesktop = com.mathworks.mde.desk.MLDesktop.getInstance;
    jDesktop.getMainFrame.setTitle([host(1:end-1),': ',numCores,' ',cpuClock,' , ',memSize,' GB RAM']);

end

