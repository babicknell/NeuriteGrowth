function Fig2BC()

load('./grad_data.mat')

cn = [0.001, 0.01, 0.03, 0.1, 0.3, 1, 3, 10, 30, 100];
cf = [0.92, 1.1, 1.5, 2.1, 1]; %correction factors from Mortimer et al. (2009)
s = [0.12 , 0.18, 0.24, 0.3, 0];

outgrowth{length(s), length(cn)} = [];
bias{length(s), length(cn)} = [];
n_val_og(length(s), length(cn)) = 0;
n_val_bs(length(s), length(cn)) = 0;

for j= 1:length(s)
    for k = 1:length(cn)
        jk = find((grad_data.NGF == cn(k)) & (grad_data.gradient == s(j)));
        outgrowth{j,k} = grad_data.averageOutgrowth(jk);
        bias{j,k} = grad_data.directionalBias(jk);
        bias{j,k}(outgrowth{j,k} < 100)=[];
        n_val_og(j,k) = length(outgrowth{j,k});
        n_val_bs(j,k) = length(bias{j,k});
    end
end

og_mean = cell2mat(cellfun(@mean, outgrowth, 'uni', 0));
og_std = cell2mat(cellfun(@std, outgrowth, 'uni', 0));
bs_mean = cell2mat(cellfun(@mean, bias, 'uni', 0));
bs_std = cell2mat(cellfun(@std, bias, 'uni', 0));


%% plotting %%
figure(1)
hold on
figure(2)
hold on
lin = {'b-','r-','g-','k-','m'};

for ind = 1:5

    figure(1)
    h(ind) = errorbar(cn.*cf(ind), og_mean(ind,:), og_std(ind,:)./sqrt(n_val_og(ind,:)),...
        lin{ind}, 'MarkerSize', 8, 'LineWidth', 2);
    
    figure(2)
    errorbar(cn.*cf(ind), bs_mean(ind,:), bs_std(ind,:)./sqrt(n_val_bs(ind,:)),...
        lin{ind},'MarkerSize',8,'LineWidth',2)
    
end


figure(1)
xlabel('NGF (nM)', 'Interpreter', 'Latex', 'FontSize', 22)
ylabel('radial outgrowth ($\mu{m}$)', 'Interpreter', 'Latex', 'FontSize', 22)
set(gca,'XScale', 'log', 'FontSize', 20)
set(gca,'xtick',[0.001,0.01,0.1,1,10,100],'xticklabel',[0.001,0.01,0.1,1,10,100]); 
legend(h,{'0.12 %','0.18 %','0.24 %','0.3 %','0 %'})
axis([0.0002,500,0,800])

figure(2)
xlabel('NGF (nM)', 'Interpreter', 'Latex', 'FontSize', 22)
ylabel('directional bias', 'Interpreter', 'Latex', 'FontSize', 22)
set(gca,'XScale', 'log', 'FontSize', 20)
set(gca,'xtick',[0.001,0.01,0.1,1,10,100],'xticklabel',[0.001,0.01,0.1,1,10,100]); 
axis([0.0002,500,-0.1,0.2])