function FigS3B()

%% 48 h %%
load('./single_cells.mat')

%indices for each concentration%
cn = [0.3, 10];
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

% plotting
figure()
hold on

h = boxplot(n_length,'positions',[1, 2],'labels',cn);
get(h,'tag');
for k = 0:1
set(h(5 + k*7),'color','r')
set(h(6 + k*7),'color','k')
end

plot([1, 2],nanmean(n_length),'ks','markersize',8,'markerfacecolor','k')

%% 96 h %%
clear all
load('./single_cells_96h.mat')

%indices for each concentration%
cn = [0.3, 10];
IND{length(cn)}=[];
for k = 1:length(cn)
        IND{k} = find(single_cells_96h.NGF == cn(k));
end
n_vals = cell2mat(cellfun(@(x) max(length(x)), IND,'uni',0));

%pad array with nans and fill with data
n_length = nan(max(n_vals),length(cn));
for k = 1:length(cn)
    n_length(1:n_vals(k),k) = single_cells_96h.length(IND{k});
end

% plotting
hold on

h = boxplot(n_length,'positions',[3, 4],'labels',cn);
get(h,'tag');
for k = 0:1
set(h(5 + k*7),'color','r')
set(h(6 + k*7),'color','k')
end

plot([3, 4],nanmean(n_length),'ks','markersize',8,'markerfacecolor','k')

axis([0,5,0,1750])
set(gca,'Xtick',[1,2,3,4],'Xticklabels',[0.3,10,0.3,10],'FontSize',22)
ylabel('length of longest neurite ($\mu m$)','Interpreter','latex','Fontsize',22)
xlabel('NGF (nM)','Interpreter','latex','Fontsize',22)
box off
