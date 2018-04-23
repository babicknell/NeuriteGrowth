function Fig1F()

load('./single_cells.mat')

%indices for each concentration%
cn = [0.003,0.01,0.03,0.1,0.3,1,10];
IND{length(cn)}=[];
for k = 1:length(cn)
        IND{k} = find(single_cells.NGF == cn(k));
end
n_vals = cell2mat(cellfun(@(x) max(length(x)), IND,'uni',0));

%pad array with nans and fill with data
n_length = nan(max(n_vals),length(cn));
for k = 1:length(cn)
    n_length(1:n_vals(k),k) = single_cells.length(IND{k});
end

%% plotting
figure()
hold on

h = boxplot(n_length,'positions',log10(cn),'labels',cn);
get(h,'tag');
for k = 0:6
set(h(5 + k*7),'color','r')
set(h(6 + k*7),'color','k')
end

plot(log10(cn),nanmean(n_length),'ks','markersize',8,'markerfacecolor','k')

axis([-3,2,0,1200])
xlabel('NGF (nM)','Interpreter','latex','Fontsize',22)
set(gca,'Xtick',[-3,-2,-1,0,1,2],'Xticklabels',[0.001,0.01,0.1,1,10,100],'FontSize',22)
ylabel('length of longest neurite ($\mu m$)','Interpreter','latex','Fontsize',22)
box off
