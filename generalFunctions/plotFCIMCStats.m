%% import FCIMCStats data
% importdata() does not work properly all time
% data = importdata('FCIMCStats');

% decide if data from NECI of kNECI
real = 1;

if real == 1
    nColumns = 31;
else
    nColumns = 28;
end

% use generic %g datatype identifier an let matlab decide
% dataFormat = repmat('%n ',1,nColumns);
% does not work with -Infinity entries: read them as strings an convert
dataFormat = [repmat('%n ',1,10),'%s ',repmat('%n ',1,11),'%s ',repmat('%n ',1,8)];

% open file
fileID = fopen('FCIMCStats');

% read first three rows only containing strings which can be used as flags
statsInfo = fgetl(fileID);
columnNames = fgetl(fileID);

% read in remaining data
data = textscan(fileID,dataFormat,'treatAsEmpty',{'Infinity','-Infinity','Infi','-Infi'},'delimiter',' ','commentStyle', '#','MultipleDelimsAsOne',1);

fclose(fileID);

% data
% time = data{1}*0.0001;
time = [1:numel(data{1})]*0.0001;
walkerPop = data{5};
HFweight = data{26};
averageShift = data{10};
projectedEnergy = data{9};

% convert energy string entries back to numericals
totEnergy = str2num(char(regexprep(data{23},'Infinity','inf')));
projEnCycle = str2num(char(regexprep(data{11},'Infinity','inf')));
% % total time somehow gets reset sometimes 
% [maxTime,indMax] = max(time);
% if indMax ~= numel(time);
%     time(indMax+1:end) = time(indMax+1:end)+maxTime;
% end

%% averaged projected energy  
%resetInd = find(data{19} <=0.001);
% energies in resets
% naiveEstimates = zeros(numel(resetInd),1);
% for i = 1:numel(resetInd)-1
%     naiveEstimates(i) = mean(totEnergy(resetInd(i):resetInd(i+1)));
% end
% naiveEstimates(end) = mean(totEnergy(resetInd(end):end));
% 
% save('naiveEstimate.mat','naiveEstimates')
% walkerPop = [100e3,500e4,1e6,2e6];


%% plots
defaultAxesPos = [0.14      0.11672        0.82     0.78];

%% number of walkers
f = figure(1);
plot(time,walkerPop,'linewidth',2)
xlabel('\bf $\tau$','fontsize',50);%ylabel('\bf walker population','fontsize',44)
% ylim([0,12e6])
% xlim([0,3.5])
set(gca,'fontsize',50,'position',defaultAxesPos)
% set(gca,'ytick',[1:6]*1e7)
saveas(f,'walkerPop.fig')
save2pdf('walkerPop',f)


%%
f = figure(2);
plot(time,averageShift,'linewidth',2)
xlabel('\bf $\tau$','fontsize',50);ylabel('\bf E$_S$','fontsize',50)
set(gca,'fontsize',50)
% ylim([-4,-2])
% xlim([0,3.5])
% set(gca,'ytick',-0.6:0.1:-0.2)
saveas(f,'shift.fig')
save2pdf('shift',f)

% %% projected energy
% figure(3)
% plot(time,data{9},'linewidth',2)
% xlabel('\bf $\tau$','fontsize',40);ylabel('\bf E$_{proj}$','fontsize',40)
% set(gca,'fontsize',36)
% screen_size = get(0, 'ScreenSize');
%% total energy
f = figure(4);
plot(time,totEnergy,'linewidth',2)
xlabel('\bf $\tau$','fontsize',50);ylabel('\bf E','fontsize',50)
% ylim([-50,-40])
set(gca,'fontsize',50,'position',defaultAxesPos)
% xlim([0,11])
saveas(f,'totEnergy.fig')
save2pdf('totEnergy',f)

%% HF weight
f = figure(5);
plot(time,HFweight,'linewidth',2)
xlabel('\bf $\tau$','fontsize',50);ylabel('\bf HF weight','fontsize',50)
% xlim([0,11])
set(gca,'fontsize',50)
set(gca,'position',defaultAxesPos)
saveas(f,'HFweight.fig')
save2pdf('HFweight',f)

% %% walker pop and shift
% f = figure(6);
% [AX,H1,H2] = plotyy(time,walkerPop,time,data{2});
% set(H1,'linewidth',4);set(H2,'linewidth',4,'color','r');
% xlabel('\bf $\tau$','fontsize',50);
% set(AX(2),'ycolor','r')
% % set(AX(2),'ylim',[-2.5,0.1])
% % set(AX(1),'ylim',[0,1.1e7],'ytick',[0:0.2:1]*1e7)
% set(get(AX(1),'Ylabel'),'String','\bf Walker population','fontsize',50)
% set(get(AX(2),'Ylabel'),'String','\bf E$_S$','fontsize',50)
% set(AX,'FontSize',50)
% 
% 
% %% walker pop and averaged energy shift
% f = figure(7);
% [AX,H1,H2] = plotyy(time,walkerPop,time,data{10});
% set(H1,'linewidth',4);set(H2,'linewidth',4,'color','r');
% xlabel('\bf $\tau$','fontsize',50);
% set(AX(2),'ycolor','r')
% % set(AX(2),'ylim',[-2.5,0.1])
% % set(AX(1),'ylim',[0,2.1e6],'ytick',[0:0.25:2]*1e6)
% set(get(AX(1),'Ylabel'),'String','\bf Walker population','fontsize',50)
% set(get(AX(2),'Ylabel'),'String','\bf $\overline{E_S}$','fontsize',50)
% set(AX,'FontSize',50)
% % set(AX,'xlim',[0,820]);
