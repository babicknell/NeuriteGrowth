function Fig1C()

load('explants.mat')

%indices for each concentration%
cn = [0.003,0.01,0.03,0.1,0.3,1,10];
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

%% plotting
figure()
hold on

for k=1:length(cn)
    plot(cn(k)*ones(1,length(outgrowth{k})),outgrowth{k},'ro','MarkerSize',5)
end

plot(cn,outgrowth_mean,'k-')

for k=1:length(cn)
    errorbar(cn(k),outgrowth_mean(k),outgrowth_sd(k)/sqrt(length(outgrowth{k})),'ks',...
    'MarkerFaceColor','k','MarkerSize',8,'LineWidth',1.5)
end

set(gca,'XScale','log','FontSize',20)
set(gca,'Xtick',[0.001,0.01,0.1,1,10,100],'Xticklabels',[0.001,0.01,0.1,1,10,100],'FontSize',22)
axis([0.001,100,0,1200])  
xlabel('NGF (nM)','Interpreter','Latex','FontSize',22)
ylabel('radial outgrowth ($\mu{m}$)','Interpreter','Latex','FontSize',22)
