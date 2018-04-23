function Fig1I()

load('./bundle_widths.mat')

%pool over patches
w_thin_p = cell2mat(w_thin.seg_widths');
w_thick_p = cell2mat(w_thick.seg_widths');
w_0p3_p = cell2mat(w_0p3.seg_widths');
w_10_p = cell2mat(w_10.seg_widths');

%cumulative distributions
wid = 0:0.1:6;
P_thin=zeros(1,length(wid));
P_thick=zeros(1,length(wid));
P_0p3=zeros(1,length(wid));
P_10=zeros(1,length(wid));

for k=1:length(wid)
    P_thin(k) = sum((w_thin_p)<wid(k))/length(w_thin_p);
    P_thick(k) = sum((w_thick_p)<wid(k))/length(w_thick_p);
    P_0p3(k) = sum((w_0p3_p)<wid(k))/length(w_0p3_p);
    P_10(k) = sum((w_10_p)<wid(k))/length(w_10_p);
end

%% plotting %%
figure()
hold on
plot(wid,P_0p3,'k','Linewidth',2)
plot(wid,P_10,'r','Linewidth',2)
plot(wid,P_thin,'k:','Linewidth',2)
plot(wid,P_thick,'k--','Linewidth',2)
hold off

axis([2,5,0,1.1])
set(gca,'FontSize',20,'ytick',0:0.2:1,'xtick',2:5)
xlabel('segment width ($\mu{m}$)','Interpreter','latex','FontSize',20)
ylabel('cumulative distribution','Interpreter','latex','FontSize',20)
legend('0.3 nM NGF','10 nM NGF','thin control','thick control','Location','southeast')

