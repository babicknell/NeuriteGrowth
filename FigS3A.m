function FigS3A()

%% 48 h %%
load('./explants.mat')
   
%indices for each concentration%
cn = [0.3,10];
IND{length(cn)}=[];
for k = 1:length(cn)
        IND{k} = find(explants.NGF == cn(k));
end

%average radial outrowth grouped by concentration%
outgrowth{length(cn)}=[];
for k = 1:length(cn)
    outgrowth{k} = explants.outgrowth(IND{k});
end

%means and sdevs%
outgrowth_mean = cell2mat(cellfun(@mean,outgrowth,'uni',0));
outgrowth_sd = cell2mat(cellfun(@std,outgrowth,'uni',0));

%plotting
figure()
hold on

for k=1:length(cn)
    plot(k*ones(1,length(outgrowth{k})),outgrowth{k},'ro','MarkerSize',5)
end

plot([1,2],outgrowth_mean,'k-')

for k=1:length(cn)
    errorbar(k,outgrowth_mean(k),outgrowth_sd(k)/sqrt(length(outgrowth{k})),'ks',...
    'MarkerFaceColor','k','MarkerSize',8,'LineWidth',1.5)
end

%% 96 h %%
clear all
load('./explants_96h.mat')

%indices for each concentration%
cn = [0.3,10];
IND{length(cn)}=[];
for k = 1:length(cn)
        IND{k} = find(explants_96h.NGF == cn(k));
end

%average radial outrowth grouped by concentration%
outgrowth{length(cn)}=[];
for k = 1:length(cn)
    outgrowth{k} = explants_96h.outgrowth(IND{k});
end

%means and sdevs%
outgrowth_mean = cell2mat(cellfun(@mean,outgrowth,'uni',0));
outgrowth_sd = cell2mat(cellfun(@std,outgrowth,'uni',0));

%plotting
hold on

for k=1:length(cn)
    plot((k+2)*ones(1,length(outgrowth{k})),outgrowth{k},'ro','MarkerSize',5)
end

plot([3,4],outgrowth_mean,'k-')

for k=1:length(cn)
    errorbar(k+2,outgrowth_mean(k),outgrowth_sd(k)/sqrt(length(outgrowth{k})),'ks',...
    'MarkerFaceColor','k','MarkerSize',8,'LineWidth',1.5)
end

hold off
axis([0,5,0,1750])
set(gca,'Xtick',[1,2,3,4],'Xticklabels',[0.3,10,0.3,10],'FontSize',22)
ylabel('radial outgrowth ($\mu m$)','Interpreter','latex','Fontsize',22)
xlabel('NGF (nM)','Interpreter','latex','Fontsize',22)

