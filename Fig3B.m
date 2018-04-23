function Fig3B()

%% model %%
%parameters
At=1000;
It=1000;
kA=0.1;
kI=5;
K0=1;
K1=100;
K2=0.75;

%steady state response
A = @(s) At*s./(s+kA);
I = @(s) It*s./(s+kI);
Gstar = @(s1,s2) K0*(A(s1)./(K1+K2*A(s1)+I(s2)));

%% data %%
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

%means and sdevs normalised by maximum mean outgrowth%
outgrowth_mean = cell2mat(cellfun(@mean,outgrowth,'uni',0));
maximum = max(outgrowth_mean);
outgrowth_mean = outgrowth_mean/maximum;
outgrowth = cellfun(@(x) x./maximum, outgrowth, 'uni', 0);
outgrowth_sd = cell2mat(cellfun(@std,outgrowth,'uni',0));


%% plotting %%
figure(1)
hold on
ss = logspace(-3,2,50);
plot(ss, Gstar(ss*1.5,ss)/max(Gstar(ss,ss)), 'r-', 'LineWidth', 2)
plot(ss, Gstar(ss,ss)/max(Gstar(ss,ss)), 'k-', 'LineWidth', 2)
plot(ss, Gstar(ss,0)/max(Gstar(ss,ss)), 'k:', 'LineWidth', 2)

for k=1:length(cn)
    errorbar(cn(k),outgrowth_mean(k),outgrowth_sd(k)/sqrt(length(outgrowth{k})),'ks',...
    'MarkerFaceColor','k','MarkerSize',10,'LineWidth',1.5)
end

xlabel('NGF (nM)', 'Interpreter', 'Latex', 'FontSize', 22)
ylabel('normalised response', 'Interpreter', 'Latex', 'FontSize', 22)
set(gca, 'XScale', 'log', 'FontSize', 20)
set(gca, 'Xtick', [0.001,0.01,0.1,1,10,100], 'Xticklabels', [0.001,0.01,0.1,1,10,100]); 
axis([0.001,100,0,1.2])